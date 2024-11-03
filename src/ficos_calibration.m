function [varargout] = ficos_calibration(varargin)
%FICOS_CALIBRATION Launch main script for calibrating the FICOS simulator. 
%   ficos_calibration() Start the script for calibrating FICOS with the default
%   arguments.
%
%   See 'help calibrationOptions' for a description of the arguments.
% 
%   Copyright (c) 2013 - 2018 Jarno Vanhatalo
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
startTime = datetime('now');
varargout = cell(nargout,1);

% Structure containing calibration parameters
% These should remain constant throughout calibration
params = struct;


% Version name is used in file paths when saving results
params.version_name = 'ficosCalibration';


% Initialize Bayes optimization state parameters
iter = 0;
deltaX = Inf;
deltaF = Inf;
nmins = 1;

continuing = 0;

%% Parsing arguments

% with inputParser
iP = inputParser;
iP.CaseSensitive = false;
% Keep unmatched arguments as metadata
iP.KeepUnmatched = true;

addParameter(iP,'version_name','',@ischar)


validLoading = @(x) ismember(x, {'loading', 'loading_Aurajoki'});
addParameter(iP, 'loading', 'loading_Aurajoki', validLoading);

addParameter(iP,'continue',false,@isfile)

valid_ma = @(x) isnumeric(x) & x>=0 & mod(x,1)==0;
addParameter(iP,'moving_average', 5,valid_ma)

dataset_intens_names = {'intens','intensiiviasemat'};
dataset_intensv2_names = {'intens_v2','intens_SeiliSum','SeiliSum'};

validDatasets = [dataset_intens_names, dataset_intensv2_names];
valid_dataset = @(x) ismember(x,validDatasets) | ...
    isfile(x);
addParameter(iP,'dataset','intens_v2',valid_dataset)

% Which parametrization (set of parameters) to identify
% light-fix: Klight, LightN2fix, LighthThresh
% orginal: KlightA, KlightFC, SumAFCmax, sedrate, springEndSedCoef
% Rmax: Klight, LightN2fix, LighthThresh, Rmax=RAmax=RFCmax
% Rmax-split: Klight, LightN2fix, LighthThresh, RAmax, RFCmax
valid_paramset = @(x) ismember(x,{'light-fix','original','Rmax',...
    'Rmax-split'});
addParameter(iP,'paramset','Rmax-split',valid_paramset);

% Option to include looser buffer after boundaries; when tightening
% boundaries keep evaluations that fall inside buffer instead. Only run new
% evaluations inside tighter boundary
addParameter(iP,'buffer',false,@islogical);

% Option to run multiple simulations in parallel during optimization loop
addParameter(iP,'parallel_bo',true,@islogical);

% Likelihood function
lik_logNames = {'censoredLog','log'};
lik_log1pNames = {'censoredLog1p','log1p'};
lik_sqrtNames = {'censoredSqrt','sqrt'};
lik_censoredNames = {'censored','censoredLin','censoredNormal'};
validLiks = [lik_logNames, lik_log1pNames, lik_sqrtNames,...
    lik_censoredNames];
valid_lik_f = @(x) ismember(x, validLiks);
addParameter(iP, 'lik_f', 'censoredLog',valid_lik_f);

addParameter(iP, 'fevalMax', 3000, @(x) isscalar(x) && isnumeric(x) ...
                                        && x>1);


valid_gpstuffroot = @(s) isempty(s) || isfolder(s);
addParameter(iP, 'gpstuffroot', valid_gpstuffroot);

addParameter(iP, 'hdf5path', fullfile('data/'), @isfolder);

% Minimum and maximum amount of cores to use when running simulator in parallel
addParameter(iP, 'cores', [2, 20], @isnumeric);

%% Parsing input 

parse(iP,varargin{:});

% min and max number of parallel processes
params.cores = iP.Results.cores;
if isscalar(params.cores)
    params.cores = repmat(params.cores, 1, 2);
end


%% Preparing hdf5 input files

params.hdf5path = iP.Results.hdf5path;

% Loading file

params.loading_file = 'loading';

if(~ismember('loading', iP.UsingDefaults))
    % use base loading.hdf5 file
    params.loading_file = iP.Results.loading;
end

% Prepare loading file and check its md5sum
loading_fullpath = fullfile(params.hdf5path, [params.loading_file, '.hdf5']);
[status] = setupLoadingFile('fullpath', loading_fullpath);
if(status)
   error('Could not prepare loading file %s.', loading_file);
end
params.loading_path = loading_fullpath;

% Prepare hd_files.hdf5
% computing moving average over time for solar radiation

moving_average = iP.Results.moving_average;
params.hd_files_ma = moving_average;
params.hd_files_path = prepareHdFiles('numDays', moving_average, ...
    'hdf5dir', params.hdf5path, 'targetDir', params.hdf5path);


%% Ensure GP stuff is within Matlab path

if(ismember('gpstuffroot', iP.UsingDefaults) || ...
        ~isfolder(iP.Results.gpstuffroot))
    if(~exist('gpstuffroot', 'var'))
        % Load gpstuffroot with startup script if not already loaded
        startup
    end
    params.gpstuffroot = gpstuffroot;
else
    gpstuffroot = iP.Results.gpstuffroot;
    params.gpstuffroot = gpstuffroot;
end

%% Initializing calibration
if(~iP.Results.continue)
% Starting calibration from the beginning
  %% Setting likelihood and simulator functions
 
paramset = iP.Results.paramset;
params.thetaNames = {};
switch(paramset)
    case 'light-fix'
        params.thetaNames = {'Klight','LightN2fix','LightThres'};
    case 'original'
        params.thetaNames = {'KlightA', 'KlightFC', 'SumAFCmax', 'sedrate', 'springEndSedCoef'};
    case 'Rmax'
        params.thetaNames = {'Klight','LightN2fix','LightThres','Rmax'};
    case 'Rmax-split'
        params.thetaNames = {'Klight','LightN2fix','LightThres',...
            'RAmax','RFCmax'};
    otherwise
        error('Uknown parametrization')
end


tt_data = timetable;
if(~isfield(params,'lik_f'))
    % Load data set used in likelihood function:
    switch iP.Results.dataset
        case {'intens','intensiiviasemat'}
            % Use all data from Seili and Utö
            tt_data = load_intens(false, 'old', 'all');
        case 'intens_v2'
            % All intensive stations, but only 1.6. - 31.10. for Seili
            tt_data = load_intens(true, 'old', 'summer');
    end

    switch iP.Results.lik_f
        case lik_logNames
            params.lik_f = @(T)likelihoodWithCensored(tt_data, T, ...
                'transform','log');
        case lik_log1pNames
            params.lik_f = @(T)likelihoodWithCensored(tt_data, T, ...
                'transform','log1p');
        case lik_sqrtNames
            params.lik_f = @(T)likelihoodWithCensored(tt_data, T, ...
                'transform','sqrt');
        case lik_censoredNames
            params.lik_f = @(T)likelihoodWithCensored(tt_data, T, ...
                'transform','none');
    end
end

end
 
% Constructing prior density function
% Define priors with containers.Map, mapping parameter name to vector
% with prior mean and variance, and initial boundaries

Priors = containers.Map('UniformValues',false);
Priors('Klight') = [log(10), log(2)/2, 0.5, 40];
Priors('LightN2fix') = [log(15), log(2)/2, 0.5, 60];
Priors('LightThres') = [log(10), log(2)/2, 0.5, 40];
Priors('Rmax') = [-2, 0.4, 0.03, 0.7];
Priors('RAmax') = [-2, 0.4, 0.03, 0.7];
Priors('RFCmax') = [-2, 0.4, 0.03, 0.7];
Priors('sedrate') = [0, log(2)/2, 0.25, 4];

params.mu_prior = zeros(1,length(params.thetaNames));
params.sigma_prior = zeros(1,length(params.thetaNames));
% Limits of the parameter values
params.limits = zeros(length(params.thetaNames),2);
for i = 1:length(params.thetaNames)
    prior_vec = Priors(params.thetaNames{i});
    params.mu_prior(i) = prior_vec(1);
    params.sigma_prior(i) = prior_vec(2);
    params.limits(i,:) = prior_vec(3:4);
end

% Log-normal prior density for parameters; parameters independent of each
% other
 params.prior_f = @(x) 0.5*sum( (((log(x)-params.mu_prior)./...
     params.sigma_prior)).^2) + sum(log(x));
 params.post_f = @(x,f)params.lik_f(f)+params.prior_f(x);
 
 % Dimensionality of the parametrization
 params.d_theta = length(params.mu_prior); % simulator parameters
 % alias used in a few places
 params.d = params.d_theta;

 % log transform the limits
 % tranforms small numbers to positive numbers first
 params.transform = @(x) log(x);
 params.retransform = @(x) exp(x);

  % Original limits remain constant, ub and lb will be shrinked during calibration
 params.limits = params.transform(params.limits);
 lb = params.limits(:,1); % lower bounds
 ub = params.limits(:,2); % upper bounds

 % Buffer bounds - wider boundary, retain evaluations here for gp fitting
 % but do not consider points outside regular boundary for exploration
 
 lb_buffer = lb;
 ub_buffer = ub;
 if(iP.Results.buffer)
     params.buffer_bounds = true;
 else
     params.buffer_bounds = false;
 end


%% Bayes optimization parameters
% Minimum amount of iterations before checking other convergence criteria,
% i.e. run at least itermin iterations at each optimization phase
params.itermin = 15; 
% After itermax iterations move to next phase.
params.fevalMax = iP.Results.fevalMax;
params.itermax = params.fevalMax;
% tolX in the same scale as original parameters, using L_inf norm, i.e.
% max(x(i,:)-x(j,:)) as the distance
params.tolX = 1e-2;
params.tolF = 1e-2;
params.tolEI = 1e-3; % Stopping condition for expected improvements
params.visualization = 0;

% No simulated annealing:
params.ps = 1;

% Amount of starting locations used when looking for expected improvement maxima
params.nstarts = 10000;
% Design size if new design needed after shrinking search space

params.n_init = 50; 
% Number of design points used for mode approximation
params.n_sf = 500;
% How many standard deviations away from mode is covered in the mode approximation
% Posterior HPD Based on asymptotic posterior normality, which assumes
% large amount of data and that posterior density is increasingly concentrated
% as amount of data increases. 
params.step_size = chi2inv(0.99, params.d)/2;

% Whether to run multiple evaluations in parallel during optimization loop
params.parallel_bo = false;
if(isfield(iP.Results,'parallel_bo'))
    params.parallel_bo = iP.Results.parallel_bo;
    if(params.parallel_bo)
        disp('Running multiple evaluations in parallel in optimization loop.')
        % Smaller number of batches (iterations) when evaluating multiple
        % parametrizations in parallel
        params.itermin = 5;
        try
            % Use undocumented MATLAB feature to get number of available cores
            % https://stackoverflow.com/questions/8311426/how-can-i-query-the-number-of-physical-cores-from-matlab
            num_cores = feature('numcores');
        catch ME
            warning(ME);
            num_cores = params.cores(2);
        end
        num_cores = min(num_cores, params.cores(2));
        params.itermax = ceil(params.fevalMax/num_cores);
    end


end



%% Setting up the GP emulator
 set_gp_emu

%% Parameters for GP MCMC approximation
params.optsls = struct();
params.optsls.nomit = 0;
params.optsls.display = 2;
params.optsls.method = 'multi';
params.optsls.wsize = 0.5;
params.optsls.plimit = 5;
params.optsls.unimodal = 1;
params.optsls.mmlimits = [-1*ones(1,params.d_theta);ones(1,params.d_theta)];
params.optsls.nsamples = 20000;
params.optsls.maxiter = 100;


%% Set simulator function

% Set paths to simulator through ficosOptions
ficosOpt = ficosOptions;
params.ficosOpt = ficosOpt;

if(~isfield(params,'sim_f'))
    % Function for evaluating a single parametrization theta and returning its
    % output
 params.sim_f = @(theta, id, wid) run_sim(theta, params.thetaNames, ...
     id, wid, hd_files=params.hd_files_path, ...
     loading=params.loading_path,ficosOpt=ficosOpt);
end

if ~isfield(params, 'simParallelFcn')
    % Function handle for evaluating a batch similar to sim_f
    params.simParallelFcn = @(x) runFICOSparallel(x, ...
        params.thetaNames, 'hd_files', params.hd_files_path, ...
        'loading', params.loading_path, ...
        'cores', params.cores, 'ficosOpt',ficosOpt);
end

%%
 % Initial design size
params.n_theta_init = 50; 

if(~ismember('version_name',iP.UsingDefaults))
    % Override version name
    params.version_name = iP.Results.version_name;
end

% Storing results and Matlab log files in the same results sub-directory
params.results_path = sprintf('results/%s_%s',...
    params.version_name,getTimeStamp());
if ~isfolder('results')
  mkdir('results')
end
mkdir(params.results_path)

% Begin calibration
try
    calibration_main
catch ME
    disp(ME)
    error_path = sprintf('%s/error_%s',params.results_path,...
        getTimeStamp());
    save(error_path);
    varargout{1} = ME;
end


if(nargout>=1)&&isempty(varargout{1})
    % Success
    varargout{1} = 0;
end
if(nargout>=2)
    varargout{2} = params;
end

end