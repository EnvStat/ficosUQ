function [params] = initializeCalibrationParams(calibOpts)
    %INITIALIZECALIBRATIONPARAMS Initialize calibration settings.
    %
    %   params = initializeCalibrationParams(calibOpts) Initializes calibration
    %   settings structure 'params' that is used within the calibration scripts.
    %
    %   The argument calibOpts is a struct array containing all the arguments to
    %   the main function for calibrating FICOS, ficos_calibration(). Use the
    %   function calibrationOptions to set the arguments.
    % 
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        calibOpts struct = calibrationOptions();
    end
    if isempty(calibOpts.version_name)
        params.version_name = 'ficosCalibration';
    else
        params.version_name = calibOpts.version_name;
    end

    if isscalar(calibOpts.cores)
        params.cores = repmat(calibOpts.cores, 1, 2);
    else
        params.cores = calibOpts.cores;
    end

    %% hdf5 settings ------
    params.hdf5path = calibOpts.hdf5path;
    params.loading_file = calibOpts.loading;
    loading_fullpath = ...
        fullfile(params.hdf5path, [params.loading_file, '.hdf5']);
    [status] = setupLoadingFile('fullpath', loading_fullpath);
    if(status)
        error('Could not prepare loading file %s.', loading_file);
    end
    params.loading_path = loading_fullpath;
    moving_average = calibOpts.moving_average;
    params.hd_files_ma = moving_average;
    params.hd_files_path = prepareHdFiles('numDays', moving_average, ...
        'hdf5dir', params.hdf5path, 'targetDir', params.hdf5path);
    params.gpstuffroot = calibOpts.gpstuffroot;

    %% define parameters by the named paramset
    switch(calibOpts.paramset)
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
            error('Uknown parametrization: %', calibOpts.paramset);
    end

    % define calibration data set and the likelihood function
    switch calibOpts.dataset
        case {'intens'}
            % Use all data from Seili and Utö
            tt_data = load_intens(false, 'old', 'all');
        case {'intens_v2'}
            % All intensive stations, but only 1.6. - 31.10. for Seili
            tt_data = load_intens(true, 'old', 'summer');
        otherwise
            error('Invalid calibration dataset: %s', calibOpts.dataset);
    end
    params.tt_data = tt_data;
    params.lik_f = constructLikFcn(tt_data, calibOpts.lik_f);

    % define prior and initial limits
    [params.prior_f, params.mu_prior, params.sigma_prior, ...
        params.limits] = constructPriorFcn(params.thetaNames);
    % define (negative log) posterior density function
    params.post_f = @(x,f) params.lik_f(f)+params.prior_f(x);

    % Dimensionality of the parametrization
    params.d_theta = length(params.mu_prior);
    % Alias for the above
    params.d = params.d_theta;

    % Transformation for the simulator parameters
    params.transform = @(x) log(x);
    params.retransform = @(x) exp(x);
    % transform the limits
    params.limits = params.transform(params.limits);

    params.buffer_bounds = calibOpts.buffer;

    %% Bayes Optization parameters

    % Minimum and maximum amount of iterations and function evaluations
    params.itermin = 15;
    params.fevalMax = calibOpts.fevalMax;
    % After itermax iterations move to next phase.
    params.itermax = params.fevalMax;

    params.tolX = 1e-2;
    params.tolF = 1e-2;
    params.tolEI = 1e-3;
    params.visualization = 0;

    params.ps = 1;

    % Amount of starting locations used when looking for expected improvement maxima
    params.nstarts = 10000;
    % Design size if new design needed after shrinking search space

    params.n_init = 50;
    % Number of design points used for mode approximation
    params.n_sf = 500;
    % Mode approximation range
    params.step_size = chi2inv(0.99, params.d)/2;

    % Whether to run multiple evaluations in parallel during optimization loop
    params.parallel_bo = calibOpts.parallel_bo;
    if params.parallel_bo
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

    params.ficosOpt = calibOpts.ficosOpt;

    % Function for evaluating a single parametrization theta and returning its
    % output
    params.sim_f = @(theta, id, wid) run_sim(theta, params.thetaNames, ...
        id, wid, hd_files=params.hd_files_path, ...
        loading=params.loading_path,ficosOpt=ficosOpt);
    if calibOpts.parallel_bo
        % Function handle for evaluating a batch similar to sim_f
        params.simParallelFcn = @(x) runFICOSparallel(x, ...
            params.thetaNames, 'hd_files', params.hd_files_path, ...
            'loading', params.loading_path, ...
            'cores', params.cores, 'ficosOpt',params.ficosOpt);
    end
    %% 
    % Initial design size
    params.n_theta_init = 50;

    % Storing results and Matlab log files in the same results sub-directory
    params.results_path = sprintf('results/%s_%s',...
        params.version_name,getTimeStamp());
    if ~isfolder('results')
        mkdir('results')
    end
    fprintf('Calibration results will be saved in %s\n', params.results_path);

end

