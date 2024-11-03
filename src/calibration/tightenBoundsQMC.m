% TIGHTENBOUNDSQMC Constrain optimisation bounds using Quasi Monte Carlo samples.
%
% Constrains optimization boundary hyper cube based on the probability that a 
% parametrisation in the area is within the Highest Prosterior Density region 
% for the parameters. 
% 
% Copyright (c) 2013 - 2018 Jarno Vanhatalo
% Copyright (c) 2017 - 2024 Karel Kaurila
%
visualizeQmc = false;


% start quasi random sequence from the beginning
quasiRng(-1);

% Include buffer boundaries
if(~exist("lb_buffer","var"))
    lb_buffer = lb;
end
if(~exist("ub_buffer","var"))
    ub_buffer = ub;
end

% here default to false to ensure old functionality
if(~isfield(params,'buffer_bounds'))
    params.buffer_bounds = false;
end
% Cutoff values for boundaries
if(~isfield(params,'boundary_cutoff'))
    params.boundary_cutoff = 10;
end
if(~isfield(params,'buffer_cutoff'))
    params.buffer_cutoff = 13;
end

%%
% Optimize GP hyperparameters and transform variables

[xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,y,...
    gp,opt,params.transform,lb,ub);

dim_x = size(xmin,2);

fh_e = @(xp) gp_cov(gp,xp,xnm)*a;
fh_g = @(xp) (grad(xp, fh_e));
fh_eg = @(xp) deal(fh_e(xp), fh_g(xp));

xstart = xnm(find(ynm==min(ynm),1),:);
ystart = min(ynm);

xmode = optimf(fh_eg, xstart, [], [], [], [], ...
    -1*ones(params.d,1), ones(params.d,1), [], optacq);
emode = fh_e(xmode);
fprintf('GP mode: g(');
fprintf('%1.3g, ', xmode(1:(end-1)) );
fprintf('%1.3g) = %1.5g\n', xmode(end), emode*stdynm+mynm);

%%

fprintf('Observed mode: h(')
fprintf('%1.3g, ', xstart(1, 1:(dim_x-1)) );
ymin = ystart*stdynm+mynm;
fprintf('%1.3g) = %1.5g\n', xstart(1, dim_x), ystart*stdynm+mynm);


if(~isfield(params, 'nQmcBounds'))
    % Amount of quasi Monte Carlo points to use for HPD region
    nFill = 500;
    params.nQmcBounds = min(nFill*3^dim_x, 1e5);
end
if(~isfield(params,'pHPD'))
    % Probabilty treshold used for high posterior density region
    params.pHPD = 0.9;
end
if(~isfield(params,'pHpdBuffer'))
    % Larger probability threshold for buffer values
    % Values within these bounds are included in GP fitting,
    % but optimization will only focus on the tigher HPD region
    params.pHpdBuffer = 0.995;
end
if(~isfield(params,'pBounds'))
    % How certain we want to be about ruling out points
    params.pBounds = 0.95;
end

%% Use GP fit to estimate probability that given point is within HPD

wQmc = quasiRng(params.nQmcBounds, size(xmode,2));

% Mapping points over slightly larger area than original so that
% 1) the tested points will include current bounds
% 2) we can test whether the previous bounds were too strict
p_extend = 0.05;
xQmc = -1*(1+p_extend)*ones(size(wQmc))+2*(1+p_extend).*wQmc;

% for a Gaussian rv, 2*lpdf difference from mode follows chi-squared
% distribution
cThresh = chi2inv(params.pHPD, dim_x)/2;
cBuffer = chi2inv(params.pHpdBuffer, dim_x)/2;
fprintf('Threshold for HPD inclusion: %1.2g\n', cThresh);
fprintf('Threshold for including in GP fit: %1.2g\n', cBuffer);

[mQmc, vQmc] = gp_pred(gp, xnm, ynm, xQmc);
sQmc = sqrt(vQmc).*stdynm;
deltaQmc = stdynm*(mQmc-min(ynm));
if(useLogOfLogPosterior)
    yLMLmin = ymin;
    yLogMin = exp(ymin);

    % Threshold for log posterior differences induces
    % threshold for absolute log minus log posterior densities
    % Not very sensitive to difference treshold c,
    % since usually log(yLogMin+c) ~= log(yLogMin)
    cLML = log(cThresh + yLogMin);
    cBufLML = log(cBuffer + yLogMin);

    pQmc = normcdf(cLML-mynm, mQmc.*stdynm, sQmc);
    pQmcBuffer = normcdf(cBufLML-mynm, mQmc.*stdynm, sQmc);
else
    pQmc = normcdf(cThresh, deltaQmc, sQmc);
    pQmcBuffer = normcdf(cBuffer, deltaQmc, sQmc);
end


% Define bounds so that parameterizations outside bounds have high
% probability of NOT being inside HPD,
% i.e. we want a minimal probability to rule out a potential region too
% early
withinBounds = (pQmc > 1-params.pBounds);
withinBuffer = (pQmcBuffer > 1-params.pBounds);

xBounds = xQmc(withinBounds,:);
xBuffer = xQmc(withinBuffer,:);

%% Estimate new bounds based on the bounding box of the QMC samples

lb_new = lb; ub_new = ub;
lb_buffer_new = lb_buffer;
ub_buffer_new = ub_buffer;

% normalized x to transformed x
xnm2trx = @(x, lb, ub) lb(:)'+(x+1)/2.*(ub(:)'-lb(:)');
% transformed to normalized ([lb, ub] -> [-1, 1])
trx2xnm = @(x, lb, ub) 2*(x-lb(:)')./(ub(:)' - lb(:)')-1;
% retransform back to original scale
fRt = @(x) params.retransform(x);
% original to transformed
fTr = @(x) params.transform(x);
% normalized to original scale
xnm2x = @(xnm, lb, ub) fRt(xnm2trx(xnm, lb, ub));

% print boundary change
printFcn = @(lb1, ub1, lb2, ub2) fprintf(...
    'Parameter %d: [%5.3g, %5.3g] -> [%5.3g, %5.3g]\n', ...
    [(1:length(lb1))' lb1(:) ub1(:) lb2(:) ub2(:)]');
printOrigFcn = @(lb1, ub1, lb2, ub2) ...
    printFcn(fRt(lb1), fRt(ub1), fRt(lb2), fRt(ub2));

if size(xBuffer,1) >= 2
    % new buffer in transformed scale
    [lb_buffer_new, ub_buffer_new] = bounds(xnm2trx(xBuffer, lb, ub));
    lb_buffer_new = lb_buffer_new(:);
    ub_buffer_new = ub_buffer_new(:);
    % print in original scale
    fprintf('Updating buffer:\n')
    printOrigFcn(lb_buffer, ub_buffer, lb_buffer_new, ub_buffer_new);
else
    error('Not enough samples inside new buffer.')
end

if size(xBounds,1) >= 2
    [lb_new, ub_new] = bounds(xnm2trx(xBounds, lb, ub));
    lb_new = lb_new(:);
    ub_new = ub_new(:);
else
    warning(['Not enough samples within constrained bounds. ' ...
        'Using buffer as bounds instead.']);
    lb_new = lb_buffer_new;
    ub_new = ub_buffer_new;
end
% Boundaries cannot expand from original
lb_new = max(lb_new, lb);
ub_new = min(ub_new, ub);

fprintf('Updating bounds:\n');
printOrigFcn(lb, ub, lb_new, ub_new);

lb = lb_new;
ub = ub_new;
lb_buffer = lb_buffer_new;
ub_buffer = ub_buffer_new;

%% Check which function evaluations remain inside buffer

iObsInBuffer = true(size(y));
lbBufRt = fRt(lb_buffer_new);
ubBufRt = fRt(ub_buffer_new);
for i = 1:dim_x
    iOut = x(:, i) > ubBufRt(i) | x(:,i) < lbBufRt(i);
    iObsInBuffer(iOut) = false;
end

x = x(iObsInBuffer,:);
y = y(iObsInBuffer);