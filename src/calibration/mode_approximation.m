%% MODE_APPROXIMATION Laplace approximation near mode.
% mode_approximation Generate a space filling design near posterior mode in
% order to obtain more information about posterior shape.
% 
% Copyright (c) 2013 - 2018 Jarno Vanhatalo 
% Copyright (c) 2017 - 2024 Karel Kaurila

%% perform new simulations near mode

if(~exist('fTimeStamp','var'))
    fTimeStamp = @() datetime('now','Format','dd-MMM-HH:mm:ss');
end

if(exist("savpath",'var'))
    % Save current state before Laplace approximation
    [saveDir, savpathStem] = fileparts(savpath);
    savNameLaplace = [savpathStem, '_pre_grid_', getTimeStamp(), '.mat'];
    savpathLaplace = fullfile(saveDir, savNameLaplace);
    save(savpathLaplace)
end
[xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,y,...
    gp,opt,params.transform,lb,ub);
fmin = min( ynm );

% Locate GP mode
xmode = gp_mode(gp, xnm, a, ynm, optacq,optimf);

% Simulate new observations near mode along principal components
fh_e = @(x) gp_pred(gp,xnm,ynm,x);

H = hessian(xmode , fh_e, []);
Sigma_orig = inv(H);
Sigma = Sigma_orig;
% Some jitter may be needed to get positive semi-definite covariance
if any(eig(Sigma_orig)<0)
    jitter = 0;
    while any(eig(Sigma)<0)
        jitter = jitter + eye(size(H,1))*0.01;
        Sigma = Sigma_orig + jitter;
    end
    warning(['ficosUQ:mode_approximation: singular Hessian. ' ...
        'Jitter of %.4f added.\n'], jitter);
end

[V,D] = eig(full(Sigma));

u_sf = quasiRng(params.n_sf, params.d);
Z_sf = 2*u_sf-ones(params.n_sf, params.d);

if(~isfield(params,'step_size'))
    params.step_size = chi2inv(0.99,size(xmode,2))/2;
end

r_grid = params.step_size/stdynm;
zgrid = repmat(xmode,params.n_sf,1) + Z_sf*(V*sqrt(D))'.*r_grid;

% Retransforming zgrid into [lb,ub]
zgrid = retransform_xy(zgrid,lb,ub,params.retransform);

% Adjust boundaries if grid goes past them
lb_preGrid = lb; ub_preGrid = ub;
lb_grid = params.transform(min(zgrid))'; 
ub_grid = params.transform(max(zgrid))';

lb = min(lb_preGrid, lb_grid);
ub = max(ub_preGrid, ub_grid);


zgrid = unique(zgrid,'rows');
params.n_sf = size(zgrid,1);

y_sf = zeros(params.n_sf,1);
y_sf_post =zeros(params.n_sf,1);

evalDesignFcn = @(x, sim_f, logOfLog) ...
    evaluate_design(x,sim_f,...
    params.post_f,tt_data,logOfLog);

if isfield(params, 'simParallelFcn')
    y_sf = evalDesignFcn(zgrid, params.simParallelFcn, useLogOfLogPosterior);
else
    y_sf = evalDesignFcn(zgrid, params.sim_f, useLogOfLogPosterior);
end

% Record where laplace grid points are in y_all
% Keep multiple rows if approximation done multiple times
if(~exist('indLaplaceGrid', 'var'))
    indLaplaceGrid = [];
end
indLaplaceGrid = [indLaplaceGrid;
                  length(y_all)+1, length(y_sf)];
% Store in x_all, y_all
x_all = [x_all; zgrid];
y_all = [y_all; y_sf];

y = [y; y_sf];
x = [x; zgrid];
