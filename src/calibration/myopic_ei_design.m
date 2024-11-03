function [x_design, varargout] = myopic_ei_design(gp, xnm, ynm,...
        bounds,options)
% MYOPIC_EI_DESIGN Construct a set of design points (myopically) by maximizing
% expected improvement.
% 
%   x_design = myopic_ei_design(gp, xnm, ynm)
%
%   Function for generating a myopic design based on expected improvement
%   Design points maximize the expected improvement conditional on previous
%   observations xnm, ynm. After first point, condition also on design points
%   x_i, y_i, i=1,...,n-1, where y_i = E[y_i| xnm, ynm]
%
%   Optional arguments:
%   
%   (positional)
%   bounds: 2 x size(xnm,2) vector with lower and upper bounds, defaults to
%   -1, 1 for each dimension if omitted
%
%   (name, value)
%   'n_design': number of design points to generate, defaults to maximum
%   parallel cores available
%
%   This strategy is a variation of the Constant Liar and Kriging Believer
%   strategies of "Kriging is Well-Suited to Parallelize Optimization" by 
%   Ginsbourger et al. (2010).
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    % Gaussian Process for the surrogate model
    gp struct
    % Parameterisations evaluated so far
    xnm {mustBeNumeric, mustBeReal}
    % Objective function values for the parameterisations xnm
    ynm (:,1) {mustBeNumeric, mustBeReal}
    % (Optional) optimization boundaries 
    % 2 x d matrix: [lower_bounds; upper_bounds]
    % default [-1, 1] for each parameter
    bounds (2,:) {mustBeNumeric, mustBeReal} = ...
        [-ones(1,size(xnm,2)); ones(1, size(xnm,2))];
    % --- (name, value) arguments ----
    % Proportion of bounds used as margin
    % Each parameter i in query points is restricted within
    % [bounds(1,i)*(1+bounds_margin), (1-bounds_margin)*bounds(2,i)] 
    options.bounds_margin (1,1) {mustBeNonnegative} = 0.02;
    % Design size, i.e. number of query points
    % default to the number of cores available or 20 if the amount of cores
    % cannot be determined
    options.n_design (1,1) {mustBePositive} = getNumCores(20);
    % how many starting points to use when optimizing expected improvement
    options.nstarts (1,1) {mustBePositive} = 1e4;
    % how many of the starting points to further refine
    % if 0, use poptimize below instead
    options.noptimize (1,1) {mustBeNonnegative} = 0;
    % alternative to noptimize, the propotion of starting points to refine
    options.poptimize (1,1) {mustBeInRange(options.poptimize, 0, 1)} = 1;
    % tolerance for inverse covariance when updating inverse covariance matrix
    % invC
    options.tolInvC (1,1) {mustBePositive} = 1e-5;
    % when covariance matrix would get close to singular, fill in remaining 
    % design points with maximin design instead of evaluating a smaller 
    % batch
    options.maximinIfSingular {mustBeNumericOrLogical} = false;
    % print messages to console
    options.verbose {mustBeNumericOrLogical} = false;
end
dim_x = size(xnm,2);

%% Options for optimizing acquisition function (expected improvement)
optimf = @fmincon;
optacqDef = optimset('fmincon');

% Tolerances for calibration are 1e-2 for both X and Fun,
% thus use slightly stricter tolerances here
optimOpt = {'GradObj','on','LargeScale','on',...
    'Algorithm','interior-point',...
    'TolFun',5e-3,'TolX',5e-3, 'Display', 'off',};
optacq = optimset(optacqDef, optimOpt{:});


%% Parse inputs

verbosePrints = options.verbose;

n_design = options.n_design;
bounds_margin = options.bounds_margin;

% bounds for starting points:
lb_start = (1-bounds_margin/2).*bounds(1,:) + bounds_margin/2.*bounds(2,:);
ub_start = bounds_margin/2.*bounds(1,:) + (1-bounds_margin/2).*bounds(2,:);


nstarts = options.nstarts;
if options.noptimize==0
    noptimize = options.poptimize * nstarts;
else
    noptimize = options.noptimize;
end
noptimize = min(noptimize, nstarts);


tolInvC = options.tolInvC;

%% Initialize
fmin = min(ynm);
xnm_ei = xnm;
ynm_ei = ynm;
x_design = zeros(n_design, dim_x);
ey_design = zeros(1,n_design);
ey_design_trunc = ey_design;
ei_design = zeros(1,n_design);

% All evaluated EIs, used for stopping condition
EI_all = zeros(n_design, nstarts);
x_cand_all = zeros(n_design, nstarts, size(xnm,2));

% Initialize invC and a that are used in calculating expected improvement
% Use noiseless covariance, since we are emulating a deterministic function
C = gp_trcov(gp, xnm_ei);
invC = inv(C);
a = C\ynm_ei;

% Define original EI function and use at the end to check for convergence
a_orig = a;
invC_orig = invC;
fh_ei_orig = @(x_new)expectedimprovement_eg(x_new, gp, xnm, a, invC_orig, fmin);
for i = 1:n_design
    if verbosePrints
        fprintf('\nQuery point %d: ', i)
    end
    fh_eg = @(x_new) expectedimprovement_eg(x_new, gp, xnm_ei, a, invC, fmin);

    % Starting points:
    x_start = quasiRng(nstarts, dim_x);

    % Map inside bounds
    x_start = lb_start+x_start.*(ub_start-lb_start);

    % Find maximum EI from starting points
    x_cand_iter = x_start;
    EIs = fh_eg(x_cand_iter);

    % Further optimize top n candidates
    [EIs, IndEIs] = sort(EIs);
    x_cand_iter = x_cand_iter(IndEIs,:);
    x_cand_optim = x_cand_iter(1:noptimize,:);
    x_cand = x_cand_optim;
    EI_all(i, :) = EIs;
    parfor(j = 1:noptimize)
        try
            x_cand(j,:) = optimf(fh_eg, x_cand_optim(j,:), [], [], [], [], ...
                bounds(1,:), bounds(2,:), [], optacq);
        catch ME
            % Retain starting point if refinement fails
            x_cand(j,:) = x_cand_optim(j,:);
        end
    end
    % Update matrix of all candidate points
    x_cand_tmp = x_cand_optim;
    if(any(isnan(x_cand)))
        % If optimization failed, use the starting point instead
        nFailStart = sum( any(isnan(x_cand),2));
        x_cand_tmp(~isnan(x_cand(:,1)),:) = x_cand(~isnan(x_cand(:,1)),:);
    end
    x_cand_iter(1:noptimize,:) = x_cand_tmp;
    x_cand_all(i, :, :) = x_cand_iter;

    % Use all points as candidates, some of which are further refined
    % through optimization, while others are the starting poiht
    x_cand = unique(x_cand_iter, 'rows', 'stable');
    EIs = fh_eg(x_cand);

    n_cand = length(EIs);
    [sortedEIs, IndEIs] = sort(EIs);
    EI_all(i, 1:n_cand) = sortedEIs;
    % Select candidate that maximizes EI
    x_design(i,:) = x_cand(IndEIs(1),:);

    % Update invC and a with the chosen candidate:
    % Inversion of a partitioned matrix,
    % (Rasmussen & Williams, 2006, p.201, equation (A.12))
    % this itself based on (Press et al., 1992, p.77)
    % here we alredy have know the inverse of the largest matrix, C
    % and only have to update one row and column
    k = gp_cov(gp, xnm_ei, x_design(i,:));
    c = gp_trvar(gp, x_design(i,:));
    invM = c - k'*invC*k;

    if(invM<tolInvC)
        % Conditional variance for candidate point close to zero, i.e.
        % the point is probably close to some earlier point
        % Try other candidates until suitable point is found
        i_cand = 2;
        found_alternative = false;
        while(invM<tolInvC && i_cand<=n_cand)
            x_design(i,:) = x_cand(IndEIs(i_cand),:);
            k = gp_cov(gp, xnm_ei, x_design(i,:));
            c = gp_trvar(gp, x_design(i,:));
            invM = c - k'*invC*k;
            i_cand = i_cand+1;
            if(invM >= tolInvC)
                found_alternative = true;
            end
        end
        if(~found_alternative)
            % All candidates would produce singular covariance matrix when
            % updating

            if options.maximinIfSingular
                nFill = n_design - i + 1;

                x_design(i:n_design, :) = maximinDesign(x_cand, nFill, xnm_ei);
            else
                % Stop looking for candidates here, return only the candidates
                % found so far
                n_design = i-1;
                x_design = x_design(1:n_design,:);
                ei_design = ei_design(1:n_design);
                ey_design = ey_design(1:n_design);
                ey_design_trunc = ey_design_trunc(1:n_design);
                x_cand_all = x_cand_all(1:n_design,:,:);
                EI_all = EI_all(1:n_design,:);
            end
            break
        end
        % Otherwise continue with the first candidate for which the updated
        % covariance is positive definite
    end
    if verbosePrints
        fprintf('[');
        fprintf('%1.3g, ', x_design(i,1:(dim_x-1)));
        fprintf('%1.3g]', x_design(i,dim_x));
    end
    M = 1/invM;
    k_inv = -invC*k*M;
    invC_tmp = invC+invC*k*M*k'*invC;

    invC = [invC_tmp, k_inv;
        k_inv', M];

    % Point estimate for y_{n+1} = E[y_{n+1}|y_n] = cov(y_{n+1},y_n)*a
    ey_design(i) = k'*a;
    % Truncate to fmin in case mean estimate is smaller than current min
    % to repel from earlier candidates. Similar to Constant Liar
    % strategies in (Ginsbourger et al., 2010).
    ey_design_trunc(i) = max(ey_design(i), fmin);
    ei_design(i) = sortedEIs(1);

    if verbosePrints
        fprintf('\nEI: %1.3g; EY: %1.3g', -ei_design(i), ey_design(i));
        if(ey_design(i)<fmin)
            fprintf(' < fmin=%1.3g, truncating to fmin.\n', fmin);
        else
            fprintf('\n');
        end
    end

    % Use current design point for evaluating next points
    xnm_ei = [xnm_ei; x_design(i,:)];
    ynm_ei = [ynm_ei; ey_design_trunc(i)];
    a = invC*ynm_ei;
end

if(nargout>=1)
    % Calculate based on the original evaluations:
    ei_design = fh_ei_orig(x_design);
    varargout{1} = ei_design;
end
if(nargout>=2)
    ey_design = gp_pred(gp, xnm, ynm, x_design);
    varargout{2} = ey_design;
end
if(nargout>=3)
    for i = 1:n_design
        x_tmp = squeeze(x_cand_all(i, :, :));
        EI_all(i, :) = fh_ei_orig(x_tmp);
    end
    varargout{3} = EI_all;
end
if(nargout>=4)
    % Return query points correponding to EI_all
    varargout{4} = x_cand_all;
end

end
