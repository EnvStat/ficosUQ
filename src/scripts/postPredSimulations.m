function [ficosOut] = postPredSimulations(options, parOpt, outputOpt)
    % POSTPREDSIMULATIONS Launch posterior predictive simulations.
    %
    % postPredSimulations('sampleFile', calibrationResults) Runs posterior
    % predictive simulations using samples from the calibration results saved in
    % the file calibrationResults (defaults to a file with example
    % results).
    % 
    % (Name, value) arguments:
    % 'sampleFile' - Path to .mat file with calibration results
    % 'nSamples' - number of posterior samples to use (default 200)
    % 'outputRoot' - where to save the results. By default, the simulation
    % results are saved in a subfolder in the same folder that contained the
    % calibration results.
    % 'cores' - minimum and maximum cores used for running the simulations in
    % parallel, default [5, 10].
    %
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.sampleFile {mustBeTextScalar} = ...
            fullfile('results/calibration_example/calibration_result.mat');
        options.nSamples (1,1) {mustBePositive, mustBeInteger} = 200;
        options.thinning (1,1) {mustBeNonnegative, mustBeInteger} = 30;
        % paths to simulator
        options.ficosOpt struct = ficosOptions;
        % hdf5 input file paths
        % if empty (default), use the paths as defined in the results
        options.hdf5opt struct = struct;
        % parallel options
        parOpt.cores {mustBePositive, mustBeInteger} = [5, 10];
        % output options --------------------
        % apply this function to each simulations output folder
        % input argument is path to the output folder
        % default: print timestamp
        outputOpt.outputFcn ...
            {mustBeA(outputOpt.outputFcn, 'function_handle')} = ...
            @(x) fprintf('%s : Simulation finished. \n', getTimeStamp());
        % simulation outputs are under this folder,
        % defaults to the be under the calibration results
        outputOpt.outputRoot {mustBeTextScalar} = '';
        % embed metadata into results
        outputOpt.embedMetadata {mustBeNumericOrLogical} = true;
        outputOpt.scenarioOpt struct = scenarioOptions();
    end

    %% Load posterior parameter samples
    if ~isempty(options.sampleFile) && isfile(options.sampleFile)
        fileVars = who('-file', options.sampleFile);
        if ~ismember('params', fileVars)
            error('variable ''params'' required to get simulator parameter names.')
        end
        load(options.sampleFile, 'params');


        if ismember('energies', fileVars)
            % load log proposal weights
            load(options.sampleFile, 'energies');
            energiesLoaded = true;
        else
            energiesLoaded = false;
        end


        if ismember('samples', fileVars)
            % load thinned samples directly
            load(options.sampleFile, 'samples');
            % 
            if size(samples, 1) > options.nSamples
                % use last n samples
                nSampTotal = size(samples,1);
                iSamp = (nSampTotal-options.nSamples+1):nSampTotal;
                samples = samples(iSamp, :);
                if energiesLoaded
                    energies = energies(iSamp);
                end
            end
        elseif ismember('ws_rt', fileVars)
            % subsample from raw samples instead
            load(options.sampleFile, 'ws_rt');
            nSampRaw = size(ws_rt,1);
            indSamp = 1:nSampRaw;
            % thinning
            if options.thinning
                indSamp = indSamp(1:options.thinning:end);
            end
            % take last n samples after thinning
            indSamp = indSamp((end-options.nSamples+1):end);
            samples = ws_rt(indSamp,:);
            if energiesLoaded
                energies = energies(indSamp);
            end
        else
            error('No posterior samples found.');
        end
    else
        error('No valid file for posterior samples given.');
    end



    %% Prepare for running simulations

    if isempty(outputOpt.outputRoot) || ~isfolder(outputOpt.outputRoot)
        outputOpt.outputRoot = fullfile(params.results_path, ...
            'postPredSimulations', getTimeStamp());
    end
    mkdir(outputOpt.outputRoot);

    if energiesLoaded
        % save sample information in a table
        createSampleTable(samples, 'energies',energies, ...
            'savePath', outputOpt.outputRoot);
    end

    % save calibration settings to a file
    save(fullfile(outputOpt.outputRoot, 'calibrationParams.mat'), 'params');

    ficosOpt = options.ficosOpt;
    hdf5Opt = options.hdf5opt;
    % defaults in case hdf5 paths not defined anywhere else
    hdf5OptDefault = hdf5options();
    % include hdf5 paths within the ficos options:
    if ~isempty(hdf5Opt) && isfield(hdf5Opt, 'hd_files')
        ficosOpt.hd_files = hdf5Opt.hd_files;
    elseif isfield(params, 'hd_files')
        ficosOpt.hd_files = params.hd_files;
    else
        ficosOpt.hd_files = hdf5OptDefault.hd_files;
    end

    if ~isempty(hdf5Opt) && isfield(hdf5Opt,'loading')
        ficosOpt.loading = hdf5Opt.loading;
    elseif isfield(params, 'loading')
        ficosOpt.loading = params.loading;
    else
        ficosOpt.loading = hdf5OptDefault.loading;
    end

    %% Run simulations

    ficosOut = runFICOSparallel(samples, params.thetaNames, ...
        cores=parOpt.cores, ...
        ficosOpt=ficosOpt, ...
        outputFcn=outputOpt.outputFcn, ...
        outputRoot=outputOpt.outputRoot);
    pause(1);
    fprintf('%s : Batch %s finished!\n', ...
        getTimeStamp(), outputOpt.outputRoot);

    %% Embed metadata in result files
    if outputOpt.embedMetadata
        fprintf('Embedding metadata in results files.\n');
        outdirs = dir(outputOpt.outputRoot);
        outdirs = outdirs(~ismember({outdirs.name}, {'.', '..'}));
        outdirs = outdirs([outdirs.isdir]);
        baseScenario = outputOpt.scenarioOpt;
        ndirs = length(outdirs);
        for i = 1:ndirs
            resultFile = fullfile(outdirs(i).folder, outdirs(i).name, ...
                'results.hdf5');
            if ~isfile(resultFile)
                warning('Could not find results file: %s', resultFile);
                continue
            end
            sampleID = str2double(outdirs(i).name);
            scenDef = baseScenario;
            scenDef.sampleID = sampleID;
            writeScenarioMetadata(resultFile, scenDef);
        end
    end
end



