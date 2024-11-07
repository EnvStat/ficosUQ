function [scenDefs] = runCatchmentScenarios(postPredOpt, options)
    %RUNCATCHMENTSCENARIOS Launch catchment load reduction scenario predictions. 
    %   runCatchmentScenarios('sampleFile', file) Run posterior predictive 
    %   simulations for catchment area load reduction scenarios using the
    %   calibraiton results in given file.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    arguments
        postPredOpt.sampleFile {mustBeTextScalar} = ...
            fullfile('results/calibration_example/calibration_results.mat');
        postPredOpt.nSamples (1,1) {mustBePositive, mustBeInteger} = 200;
        postPredOpt.thinning (1,1) {mustBeNonnegative, mustBeInteger} = 30;
        % paths to simulator
        postPredOpt.ficosOpt struct = ficosOptions;
        % parallel options
        postPredOpt.cores {mustBePositive, mustBeInteger} = [5, 10];
        % relative catchment area loadings in each scenario
        options.catchmentLoads (1,:) {mustBeInRange(options.catchmentLoads,...
            0,1)} = (1:-0.2:0.2);
        % where to save the results 
        options.scOutRoot {mustBeTextScalar} = '';
    end

    % scenario definitions
    cLoads = options.catchmentLoads;
    nScen = length(cLoads);
    for i = 1:nScen
        scenDefs(i) = scenarioOptions('scenarioID', i, ...
            'catchmentLoad',cLoads(i));
    end


    if isempty(options.scOutRoot)
        calibResDir = fileparts(postPredOpt.sampleFile);
        options.scOutRoot = fullfile(calibResDir, 'catchmentScenarios', ...
            getTimeStamp());
    end
    mkdir(options.scOutRoot);
    fprintf('Storing scenario results in %s.\n', options.scOutRoot);

    %% Create scenario loading files
    load(postPredOpt.sampleFile, 'params');
    loadingSource = params.loading;
    scenLoadingPath = fullfile(options.scOutRoot, 'loading_scenario');
    loadingFiles = compose('%s_%d.hdf5', scenLoadingPath, 1:nScen);
    for i = 1:nScen
        change_loadings(loadingSource, 'copy', loadingFiles{i},...
            'int', scenDefs(i).internalLoad, ...
            'catchment', scenDefs(i).catchmentLoad, ...
            'atm', scenDefs(i).atmosphericLoad, ...
            'point', scenDefs(i).pointLoad);
    end
    hdf5Opt = struct();
    hdf5Opt.hd_files = params.hd_files;
    
    %% Run posterior predictions
    for i = 1:nScen
        scenName = sprintf('scenario_%d', i);
        scenDir = fullfile(options.scOutRoot, scenName);
        mkdir(scenDir);
        hdf5Opt.loading = loadingFiles{i};
        postPredArgs = namedargs2cell(postPredOpt);
        postPredSimulations(postPredArgs{:}, 'hdf5opt', hdf5Opt, ...
            'outputRoot', scenDir, 'embedMetadata', true, ...
            'scenarioOpt', scenDefs(i));
    end
end

