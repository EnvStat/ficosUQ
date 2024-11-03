
%% Including external Matlab code ---------------
% Installed as submodules

% gpstuff 
gpstuffroot = fullfile('src/submodules/gpstuff/');

addpath(fullfile(gpstuffroot, 'diag'));
addpath(fullfile(gpstuffroot, 'dist'));
addpath(fullfile(gpstuffroot, 'gp'));
addpath(fullfile(gpstuffroot, 'misc'));
addpath(fullfile(gpstuffroot, 'optim'));
addpath(fullfile(gpstuffroot, 'tests'));

% External functions for approximation gradient and hessian with finite differences
% See License.txt files include in their folders
addpath src/util/external
addpath src/util/external/grad

% Color Brewer
% This product includes color
% specifications and designs developed by Cynthia Brewer
% (http://colorbrewer.org/).
% These are used for better palettes in some plots
addpath src/submodules/BrewerMap


%% Main source code
addpath src
addpath src/calibration
addpath src/scenarios
addpath src/simulator
addpath src/util
addpath src/util/options
addpath src/plotting
addpath src/posteriorPredictive

% Driver and other scripts for analyzing specific data sets or results 
addpath src/scripts




