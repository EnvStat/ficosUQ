%% OPTIMIZATION_LOOP Locate posterior mode with Bayes optimization.
% optimization_loop Script for finding posterior mode using a GP emulator and
% Bayesian optimization
%
% Copyright (c) 2013 - 2018 Jarno Vanhatalo
% Copyright (c) 2017 - 2024 Karel Kaurila

% Functions for calculating distance criteria for convergence
if ~isfield(params, 'deltaXfcn')
    params.deltaXfcn = @(x) norm(x, Inf);
end
if ~isfield(params, 'deltaXnmFcn')
    params.deltaXnmFcn = @(xnm) norm(xnm, 2);
end


if isfield(params, 'simParallelFcn')
    evalDesignFcn = @(x_new, logOfLog) evaluate_design(x_new, ...
        params.simParallelFcn, ...
        params.post_f, tt_data, logOfLog);
else
    evalDesignFcn = @(x_new) evaluate_design(x_new, ...
        params.sim_f, ...
        params.post_f, tt_data, logOfLog);
end

% whether to evaluate multiple parametrizations in parallel at each
% iteration
if(~isfield(params,'parallel_bo'))
    params.parallel_bo = false;
end

if(ind_start > numel(y))
    ind_start = numel(y);
end

if(~isfield(params, 'tolEI'))
    % Stopping condition for Expected Improvement
    params.tolEI = 1e-3;
end
if(~exist('maxEI', 'var'))
    maxEI = Inf;
end

if(~isfield(params,'tolXnm'))
    % tolerance for transformed and normalized parameter values,
    % which are used for fitting the GP
    params.tolXnm = 1e-2;
end
if(~exist('deltaXnm', 'var'))
    deltaXnm = Inf;
end


if(~exist('nBatch','var'))
    % Record batch start and end indices in y_all
    % First set is initial design to ensure correct size in the array
    nBatch = [0, length(y_all)];
end
if(~exist('indOptStart','var'))
% Record indices of batches that aries outside optimization loop
    indOptStart = length(y_all)+1;
    indTightenBounds = [];
    indLaplaceGrid = [];
end

if(~exist('fTimeStamp','var'))
    fTimeStamp = @() datetime('now','Format','dd-MMM-HH:mm:ss');
end

% Store best point from current batch for stopping condition
% Only compare the best value from each batch in case points cluster close
% together in a single batch
[ymin_prev, imin_prev] = min(y);
xmin_prev = x(imin_prev,:);

% Loop for at least itermin iterations
% After that loop until itermax iterations or until two smallest
% (log) negative log posterior densities and their parameterization are
% within tolerances (both tolerances need to be met for convergence)

optimizationDone = false;

% Convergence criteria in the while loop:
%  run at least itermin (default 5) iterations
%  and at most itermax (default 150) iterations
%  within these bounds calibration is considered converged if
%  - maximum expected improvement <= params.tolEI (default 0.001)
%  OR
%  - the two best modes are close enough to each other:
%    - values at the modes are within params.tolF
%    - parameters are within tolX in absolute scale OR
%      within tolXnm in normalized scale
while(iter <= params.itermin || ...
        (iter <= params.itermax && ...
        maxEI > params.tolEI && ...
        ~(deltaF < params.tolF && ...
         (deltaX < params.tolX || deltaXnm < params.tolXnm))))
    
    [xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,...
        y,gp,opt,params.transform,lb,ub);
    
    [fmin, iMin] = min(ynm);
    xmin = x(iMin,:);
    fh_eg = @(x_new) expectedimprovement_eg(x_new,gp,xnm,a,invC,fmin);

    if(~params.parallel_bo)
        % Run only one evaluation per iteration
        if size(xnm,1)<1
            xstart = repmat(-1*ones(1,params.d),params.nstarts,1) + ...
                repmat(2*ones(1,params.d),params.nstarts,1).*rand(params.nstarts, params.d);
        else
            xstart = repmat(-1*ones(1,params.d),2*ceil(params.nstarts/3),1) + ...
                repmat(2*ones(1,params.d),2*ceil(params.nstarts/3),1).*rand(2*ceil(params.nstarts/3), params.d);
            [~,Ibest]=sort(ynm);
            xstart3 = xnm(Ibest(1:min(floor(params.nstarts/3),length(ynm))),:);
            xstart = [xstart ; xstart3];
        end
        % Optimize the acquisition function
        x_new = [];
        if ~isempty(y)
            for s1=1:size(xstart,1)
                try
                    x_new(s1,:) = optimf(fh_eg, xstart(s1,:), [], [], [], [], ...
                        -1*ones(1,params.d), ones(1,params.d), [], optacq);
                catch err
                    x_new(s1,:) = nan(1,params.d);
                end
            end
            x_new = x_new( ~isnan(x_new(:,1)), :);
            if isempty(x_new)
                x_new = xstart(1,:);
            end
            EIs = fh_eg(x_new);
            [~, IndEIs] = sort(EIs);
            maxEI = max(-stdynm*EIs);
            x_new = x_new( IndEIs(1) , : );
        else
            x_new = xstart;
        end
        % Retransform x_new back into [lb,ub]
        x_new = retransform_xy(x_new,lb,ub,params.retransform);
       
        %calculate the objective function at new query point
        yt=[];
        yt_post=[];
        for i1 = 1:size(x_new,1)
    
            tic
            f = params.sim_f(x_new(i1,:),iter,[]);
            toc
            yt(i1) = params.post_f(x_new(i1,:),f);
        end
    else
        % Evaluate multiple points in parallel
        [xnm_new, EI_new, ey_new, EI_all] = myopic_ei_design(gp, xnm, ynm);
        x_new = retransform_xy(xnm_new, lb, ub, params.retransform);

        maxEI = max(stdynm*abs(EI_all),[],"all");
        if(useLogOfLogPosterior)
            fprintf('Expected log minus log densities:\n')
        else
            fprintf('Expected negative log densities:\n')
        end
        disp(ey_new*stdynm+mynm)

        disp('Current minimum:')
        disp(ymin)
        disp(xmin)
        yt = evalDesignFcn(x_new, false);
    end

    % Store new evaluation(s) in x_all, y_all
    nBatch = [length(y_all)+1, length(yt)];
    x_all = [x_all; x_new];
    y_all = [y_all; yt];
    if useLogOfLogPosterior
        yt = log(yt);
    end
    disp('Observed densities:')
    disp(yt(:)')
    if( any(yt < ymin))
    % new minimum
       nmins = nmins+1;
       disp('Found new minimum:')
       ymin = min(yt)
       xmin = x_new(find(yt==min(yt),1),:)
    end
    
    % New sample points
    x = [x ; x_new];    
    y = [y ; yt];

   
    % update GP fit
    [xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,...
        y,gp,opt,params.transform,lb,ub);

    % check the difference between the two best parameter combinations to
    % inform optimization stopping
    [ys, Is] = sort(y, 'ascend'); % smallest values first
    xs = x(Is,:);
    if(numel(ys)>1)

        ymin_batch = min(yt);
        xmin_batch = x_new(find(yt==ymin_batch,1),:);

        deltaF = ys(2)-ys(1);
        deltaX = params.deltaXfcn( xs(2,:) - xs(1,:));

        % Euclidian distance for transformed and normalized parameter,
        % tolerance based on GP fit stability as nearby values can cause
        % covariance to be singular
        xnms = xnm(Is,:);
        deltaXnm = params.deltaXnmFcn( xnms(2,:) - xnms(1,:));
    else
        disp('Only including parameter combinations evaluated in this optimization step in calculating deltas')
        deltaF = Inf;
        deltaX = Inf;
        deltaXnm = Inf;
    end

    % print progress reports
    fprintf('iter %d;  LPDdelta: %.5g; Xdelta: %.5g, XNMdelta: %.5g, maxEI: %.5g \n', ...
        iter, deltaF, deltaX, deltaXnm, maxEI)
    
    iter = iter+1;
    if( iter > params.itermin && ...
                    (iter>params.itermax || ...
    (deltaF < params.tolF && ...
            (deltaX < params.tolX || deltaXnm < params.tolXnm)) || ...
                    maxEI < params.tolEI ) )
       %% Save optimization state before constricting optimization boundaries
       strTimeBounds = getTimeStamp();
       savpath_bounds = sprintf('%s/bo_loop_pre_tighten_bounds_%s.mat',...
           params.results_path, strTimeBounds);
       save(savpath_bounds);
       
       save(savpath);

       %% Tighten optimization bounds
       if(~exist("indTightenBounds",'var'))
           indTightenBounds = [];
       end
       indTightenBounds = [indTightenBounds; length(y_all)];
       tightenBoundsQMC

       if (numel(y) < 50)
            disp('Less than 50 points inside new bounds => new initial design');
            if(~isfield(params,'n_init'))
                params.n_init = 50;
            end
           
            u_sf = quasiRng(params.n_init, size(x,2));
            x_sf = u_sf.*(ub-lb)'+lb';
            x_sf = params.retransform(x_sf);
            y_sf = evalDesignFcn(x_sf, useLogOfLogPosterior);
            y = [y; y_sf];
            % Store in x_all, y_all
            x_all = [x_all; x_sf];
            if(useLogOfLogPosterior)
                y_all = [y_all; exp(y_sf)];
            else
                y_all = [y_all; y_sf];
            end

            x = [x; x_sf];
       end

       %%

       if(useLogOfLogPosterior)
        if(numel(y)>0)
            % change target to minus log posterior
            y = exp(y);                     
            ymin = exp(ymin);
            ymin_batch = exp(ymin_batch);
            ymin_prev = exp(ymin_prev);
            
            deltaF = (exp(deltaF)-1)*ymin;

            [xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator( ...
                x,y,gp,opt,params.transform,lb,ub);
            fmin = min( ynm );
        end
        useLogOfLogPosterior = false;
       else
           % Optimization done
           % Only need to constrain boundaries and verify that new bounds
           % have enough evaluations
           optimizationDone = true;
       end
       
       if(optimizationDone)
           break;
       end
       
       iter = 0;
       ind_start = numel(y)+1;
       [xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(...
           x,y,gp,opt,params.transform,lb,ub);
       maxEI = Inf;
    end

end
