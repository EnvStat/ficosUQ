%% CALIBRATION_MAIN Main script for calibrating the FICOS simulator.
% calibration_main Runs scripts for each step of the calibration process and 
% saves results between steps.
% 
% The steps of the procress are:
% - initial design (initial_design.m) : Initial design before Bayes
% optimization
% - Bayes optimization (optimization_loop.m)
% - Laplace approximation near mode (mode_approximation.m)
% - MCMC sampling from GP (gp_mcmc_approximation.m)
%
% Copyright (c) 2013 - 2017 Jarno Vanhatalo
% Copyright (c) 2018 - 2024 Karel Kaurila

% The first file below will be updated between steps
savpath = sprintf('%s/calibration_results.mat',params.results_path);
save(savpath);
% Separate save for initial state (not updated after this)
save(sprintf('%s/initial_state.mat', params.results_path));

useLogOfLogPosterior=1;

% Record matlab commands and their console outputs to log files
diary off
diarypath = sprintf('%s/init_design_log.txt',params.results_path);
diary(diarypath);
diary on

% keep all evaluated densities in (x_all, y_all), which will not be
% filtered when tightening optimization boundaries
x_all = [];
y_all = [];

%% Initial design
initial_design

% Record batch sizes
nBatch = [0, length(y_all)];
% Record indices of batches that aries outside optimization loop
indOptStart = length(y_all)+1;
indTightenBounds = [];
indLaplaceGrid = [];


x_all = [x_all; x];
if(useLogOfLogPosterior)
    y_all = [y_all; exp(y)];
else
    y_all = [y_all; y];
end

[miny, indy] = min(y);
[xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,y,gp,opt,...
                                                params.transform,lb,ub);

save(savpath);
diary off

%% Bayes optimisation loop
diarypath = sprintf('%s/calibration_loop_log.txt',params.results_path);
diary(diarypath);
diary on

ymin = min(y)
% Index of the first evaluation in the Bayes optimzation process
ind_start = numel(y)+1;  
save(savpath);

piter = 1;
optimization_loop

[xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,y,gp,opt,...
    params.transform,lb,ub);

save(savpath);
diary off
diarypath = sprintf('%s/mode_approximation_log.txt',params.results_path);
diary(diarypath);
diary on

%% Laplace approximation near mode

mode_approximation
[xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(x,y,gp,opt,...
    params.transform, lb, ub);


fmin = min( ynm );
save(savpath);

diary off
diarypath = sprintf('%s/gp_mcmc_log.txt',params.results_path);
diary(diarypath);
diary on

%% GP MCMC
% Sample from emulator posterior using Markov Chain Monte Carlo

gp_mcmc_approximation

save(savpath);
diary off

fprintf('Calibration done!.\n');
fprintf('See results in %s.\n', savpath);