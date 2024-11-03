function [ficosOut] = runFICOSparallel(Theta, thetaNames, options, ...
        fixedSimOpt, outputOpt, parallelOpt, overrideOpt)
    %RUNFICOSPARALLEL Run multiple FICOS simulations in parallel.
    %   ficosOut = runFICOSparallel(Theta, thetaNames) Launches multiple FICOS
    %   simulations in parallel using Matlab's parfor, and returns the
    %   predictions in a cell array.
    %   
    %   Arguments:
    %     
    %   Theta      : (required) n x d matrix, where d is the dimension of a single
    %                 parametrisation.
    %   thetaNames : (optional) names of the parameters in a cell array. If not
    %                given, defaults to {'Klight', 'LightN2fix', ...
    %                'LightThres', 'RAmax', 'RFCmax'}.
    %   
    %  Optional named arguments
    %   'outputRoot' : folder where all outputs from this batch are stored.
    %   'cores'      : amount of parallel processes to use
    %                 (see 'NumWorkers' argument for Matlab's 'parpool').
    %   'ficosOpt'   : set paths for the FICOS simulator, see 'ficosOptions'.
    %   
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        Theta (:,:) {mustBeNumeric}
        thetaNames (1,:) cell = ...
            {'Klight', 'LightN2fix', 'LightThres', 'RAmax', 'RFCmax'};
        options.hd_files  (1,:) {mustBeTextScalar} = 'hd_files.hdf5';
        options.loading   (1,:) {mustBeTextScalar} = 'loading.hdf5';
        options.local_dir (1,:) {mustBeTextScalar} = '.';
        fixedSimOpt.simPath   (1,:) {mustBeTextScalar} = '.';
        fixedSimOpt.debugPrints (1,1) {mustBeNumericOrLogical} = false;
        fixedSimOpt.settingsSource (1,:) {mustBeTextScalar} = 'ficossettings.ini';
        fixedSimOpt.iniFile (1,:) {mustBeTextScalar} = 'wqficos.ini';
        outputOpt.outputRoot (1,:) {mustBeTextScalar} = '';
        outputOpt.rmOldOutput (1,1) {mustBeNumericOrLogical} = false;
        outputOpt.outputFcn {mustBeA(outputOpt.outputFcn, 'function_handle')} = ...
            @(x) readOutput(x);
        % Options specific to parallelization
        parallelOpt.cores (:,1) {mustBeNonnegative, mustBeInteger} = [2, 20];
        parallelOpt.parpoolOpt (:,1) cell = {}; % additional parpool options
        % override options above through struct array fields
        overrideOpt.ficosOpt struct = struct;
    end

    %% Additional input handling

    if ~isempty(overrideOpt.ficosOpt)
        % override simulator options
        fixedSimOpt = assignFields(fixedSimOpt, overrideOpt.ficosOpt);
        options = assignFields(options, overrideOpt.ficosOpt);
    end
    % Stricter input validation after assigning overriden settings
    mustBeFolder(fixedSimOpt.simPath);
    mustBeFile(fixedSimOpt.settingsSource);
    mustBeFile(fixedSimOpt.iniFile);

    if isempty(outputOpt.outputRoot)
        % Default output root folder
        outDirRoot = sprintf('batch-%s', getTimeStamp());
        outputOpt.outputRoot = fullfile('output', outDirRoot);
    end


    %% Initialize parallel pool
    cores = parallelOpt.cores;
    if ( length(cores) >= 2 )
        cores = cores(1:2);
        if ( cores(2) < cores(1) )
            error('Decreasing range: [%d %d]', cores(1), cores(2));
        end
    end
    
    % Try using undocumented Matlab feature to check for available cores
    try
        numCores = feature('numcores');
        cores = min(cores, numCores);
    catch ME
       warning('Cannot use feature numcores.')
    end

    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj = parpool(cores);
    elseif all(poolobj.NumWorkers < cores)
        % Use existing parpool, unless it has less than requested amount of
        % cores
        delete(poolobj);
        poolobj = parpool(cores);
    end

    % attach files
    filesToAttach = {fullfile(fixedSimOpt.simPath, 'wqficos'), ...
                     fixedSimOpt.iniFile, options.hd_files, ...
                     options.loading};
    
    % Fetch hdf5 input file names
    [~, hdf5Names, hdf5Ext] = fileparts({options.hd_files, options.loading});
    if ~all(strcmp(hdf5Ext, '.hdf5'))
        error('Invalid hdf5 input files.');
    end
    hd_file_path = options.hd_files;
    loading_path = options.loading;
    hd_file_name = [hdf5Names{1}, hdf5Ext{1}];
    loading_name = [hdf5Names{2}, hdf5Ext{2}];

    addAttachedFiles(poolobj, filesToAttach);
    % suppress warning
    warning('off', 'parallel:lang:pool:IgnoringAlreadyAttachedFiles');

    nBatch = size(Theta, 1);
    ficosOut = cell(nBatch, 1);

    % Other simulation options that are the same for all outputs
    fixedSimOpt = struct2opt(fixedSimOpt);
    % Suppress mkdir warnings
    warning('off', 'MATLAB:MKDIR:DirectoryExists');

    outputRoot = outputOpt.outputRoot;
    outputFcn = outputOpt.outputFcn;
    %% Run simulations
    parfor i = 1:nBatch

        % Get local input files         
        hd_files_dir = getAttachedFilesFolder(hd_file_path);
        loading_dir = getAttachedFilesFolder(loading_path);
        local_hd_files = fullfile(hd_files_dir, hd_file_name);
        local_loading = fullfile(loading_dir, loading_name);

        % Local output dir
        outputDir = fullfile(outputRoot, num2str(i));
        mkdir(outputDir);

        ficosOut{i} = run_sim(Theta(i,:), thetaNames, i, ...
            'hd_files',local_hd_files, 'loading',local_loading, ...
            'outputPath', outputDir, 'outputFcn',outputFcn, ...
            fixedSimOpt{:});
    end

end