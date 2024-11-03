% INITIAL_DESIGN Initial design before Bayes optimisation.
%
% Script for generating the initial design points and evaluating the simulator
% output and the corresponding (log) negative log posterior density for each
% design point
%
% Copyright (c) 2017 - 2024 Karel Kaurila
%

% reset sobol sequence to start from the beginning
% analoguous to resetting pseudo random number generator to default seed
[~] = quasiRng(1, params.d, 'Skip',0);

% Transforming sobol set from [0,1] to [lower bound, upper bound]
zInit = quasiRng(params.n_theta_init, params.d, 'Skip',0);
theta_init = params.limits(:,1)'+zInit.*...
    (params.limits(:,2)-params.limits(:,1))';

% including prior mean in the initial design
theta_init(end,:) = params.mu_prior;

x = params.retransform(theta_init);

% run simulator at initial design points
ymin = Inf;
fvals_min = NaN;

if ~isfield(params, 'simParallelFcn')
    y = evaluate_design(x, params.sim_f, params.post_f, tt_data, ...
        useLogOfLogPosterior, 'mode', 'legacy');
else
    y = evaluate_design(x, params.simParallelFcn, params.post_f, tt_data, ...
        useLogOfLogPosterior, 'mode', 'new');
end
