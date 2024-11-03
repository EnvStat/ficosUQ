function [varargout] = ficosCalibration(calibrationOpt)
    %FICOSCALIBRATION Calibrate parameters for FICOS simulator.
    %   ficosCalibration() Calibrate FICOS simulator using the default options.
    %   The results, including posterior parameters samples, will be saved to
    %   the to a folder under ficosUQ/results.
    %
    %   Use the function 'calibrationOptions' to define the options. Many of
    %   the options can be tab-completed, such as the path to the hdf5 input
    %   files defined by the option "hdf5path".
    %
    %   Copyright 2013 - 2018 Jarno Vanhatalo
    %   Copyright 2017 - 2024 Karel Kaurila
    %
    %   This will replace FICOS_CALIBRATION in the future, but has not been
    %   properly tested yet. For now, use FICOS_CALIBRATION instead.
    %
    arguments
        calibrationOpt struct = calibrationOptions();
    end

    %% Initialize internal settings structure 
    % used within the calibration scripts
    params = initializeCalibrationParams(calibrationOpt);

    %% Initialize remaining local variables
    
    % Initialize Bayes optimization state parameters
    iter = 0;
    deltaX = Inf;
    deltaF = Inf;
    nmins = 1;

    gpstuffroot = params.gpstuffroot;
    tt_data = params.tt_data;
    
    % optimization bounds
    lb = params.limits(:,1); % lower bounds
    ub = params.limits(:,2); % upper bounds
    % buffers (optional)
    lb_buffer = lb;
    ub_buffer = ub;
    if params.parallel_bo
        num_cores = params.cores(2);
    end
    
    % Initialize GP emulator
    set_gp_emu

    % create folder for results
    if ~isfolder(params.results_path)
        mkdir(params.results_path)
    end

    %% Begin calibration

    % wrap in try catch in case calibration is interrupted etc.
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

