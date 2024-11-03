function [cOpt, params] = calibrationOptions(calibOpts)
    %CALIBRATIONOPTIONS Define options for the main calibration function for
    %FICOS, 'ficosCalibration'. 
    %   This function is used to conveniently construct the default arguments
    %   for the ficosCalibration function and to change particular arguments.
    %   The resulting structure can be converted to the arguments used in the
    %   function call with namedargs2cell. 
    %   >> cOpt = calibrationOptions();
    %   >> cArgs = namedargs2cell(cOpt);
    %   >> ficosCalibration(cArgs{:});
    %
    %  see also FICOSCALIBRATION INITIALIZECALIBRATIONPARAMS
    %
    arguments
        % name for this calibration run, used as a name for the results 
        % folder
        calibOpts.version_name char = '';
        % which nutrient loading file to use
        % 'loading' - original loading file and
        % 'loading_Aurajoki' - loading file specific to the Aurajoki basin 
        calibOpts.loading {mustBeTextScalar, ...
            mustBeMember(calibOpts.loading, {'loading', 'loading_Aurajoki'})} = ...
            'loading_Aurajoki';
        % whether to smooth the solar radiation values in hd_files.hdf5
        % by computing a moving average
        % 0 - use the original file
        % n - (positive integer) compute a moving average over +-n days
        calibOpts.moving_average (1,1) {mustBeNonnegative, mustBeInteger} = 5;
        % which calibration data set to use
        % 'intens' - all observations from intensive stations Seili and Utö
        % 'intens_v2' - all observations from intensive stations Utö and Brändö,
        %      from Seili station use only summer observations
        calibOpts.dataset {mustBeTextScalar, ...
            mustBeMember(calibOpts.dataset, {'intens', 'intens_v2'})} = ...
            'intens_v2';
        % which set of FICOS parametes to calibrate; other parameters are
        % fixed to their default values
        % light-fix - Klight, LightN2fix, LighthThresh
        % orginal - KlightA, KlightFC, SumAFCmax, sedrate, springEndSedCoef
        % Rmax - Klight, LightN2fix, LighthThresh, Rmax=RAmax=RFCmax
        % Rmax-split - Klight, LightN2fix, LighthThresh, RAmax, RFCmax
        calibOpts.paramset {mustBeTextScalar, mustBeMember(calibOpts.paramset, ...
            {'light-fix', 'original', 'Rmax', 'Rmax-split'})} = 'Rmax-split';
        % Include a second, looser boundaries as a buffer when constraining
        % the optimization space. When using this, GP is fitted with
        % parameterisations falling inside the buffer. Calibration is always
        % constrained to the tighter boundary.
        calibOpts.buffer (1,1) {mustBeNumericOrLogical} = false;
        % Run multiple simulations in parallel during optimization loop
        calibOpts.parallel_bo (1,1) {mustBeNumericOrLogical} = true;
        % Which likelihood function to use
        calibOpts.lik_f {mustBeTextScalar, mustBeMember(calibOpts.lik_f, ...
            {'tt_v2','t', 'censoredLog','log', 'censoredSqrt','sqrt', ...
            'censored', 'censoredLin', 'censoredNormal'})} = 'censoredLog';
        % Maximum number of function evaluations to use during optmization loop
        calibOpts.fevalMax (1,1) {mustBePositive, mustBeInteger} = 3000;
        % Path to gpstuff installation
        calibOpts.gpstuffroot {mustBeTextScalar} = ...
            fullfile('src/submodules/gpstuff/');
        % Path to hdf5 input files for FICOS
        calibOpts.hdf5path {mustBeFolder} = 'data/';
        % Minimum and maximum amount of cores to use when running simulator in
        % parallel
        calibOpts.cores (1,:) {mustBePositive, mustBeInteger} = [2, 20];
        % options for the ficos simulator
        calibOpts.ficosOpt struct = ficosOptions();
    end
    % 
    cOpt = calibOpts;
    params = struct;
    if nargout >=2
        % return also the params structure that is initialized in the beginning
        % of ficos_calibration
        params = initializeCalibrationParams(cOpt);
    end
end