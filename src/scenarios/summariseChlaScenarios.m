function [tbChlaSum] = summariseChlaScenarios(resultDir, options)
    %SUMMARISECHLASCENARIOS Collect and summarise chlorophyll scenario
    %predictions.
    %   tbChlaSum = summariseChlaScenarios(resultDir) Collects and summarises
    %   chlorophyll scenario predctions in folder resultDir.
    %
    %   Optional named arguments:
    %   'weights'       : (true/false) use importance weights
    %   'savePredTable' : (true/false) save table of collected predictions to a
    %   file
    %   'saveSumTable'  : (true/false) save summary tbChlaSum to a file
    %   'plot'          : (true/false) plot distributions for mean chla based on
    %                     summary.
    %   'savePlot'      : (true/false) save the above distribution plot to file
    %   'weightBy'      : weight water formations by 'area' or by 'volume'
    %   'calibrationResults' : calibration result file to use for computing
    %                         importance weights (default:
    %                         '../calibration_results.mat')
    %   
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        resultDir {mustBeFolder}
        options.weights {mustBeNumericOrLogical} = false;
        options.savePredTable {mustBeNumericOrLogical} = true;
        options.saveSumTable {mustBeNumericOrLogical} = false;
        options.plot {mustBeNumericOrLogical} = false;
        options.savePlot {mustBeNumericOrLogical} = false;
        options.weightBy {mustBeTextScalar, mustBeMember(options.weightBy, ...
            {'area', 'volume'})} = 'area';
    end

    %% Parse inputs
    
    %% Collect all predictions
    collectOpt = struct;
    collectOpt.mode = 'hdf5';
    % Use water formations corresponding to Aurajoki estuary
    collectOpt.filtBlocks = filterBlocksArg("Aurajoki");
    % Load variables used for calibration and scenarios
    collectOpt.vars = getSimOutNames('extScenario');
    % Recurse into sub folders (scenario_1, scenario_2, ...)
    collectOpt.recurse = true;
    % Format prediction output
    collectOpt.postFcn = @(x) cleanSimOutput(x, 'new');
    collectArgs = namedargs2cell(collectOpt);
    tbPredFull = collectPostPredResults(resultDir, collectArgs{:});
    
    %% Formatting table
    oldVarNames = {'sampleID'};
    newVarNames = {'SampleID'};

    sumVars = {'sampleID', 'scenarioID', ...
        'Catchment area load', 'Point load', ...
        'Internal load', 'Atmospheric load', 'wf_id', 'chla'};
    tbPred = tbPredFull(:, sumVars);

    %% Save collected table
    if options.savePredTable
        predTbPath = getNewSaveFile(...
            fullfile(resultDir, 'scenarioPredictionTable.mat'));
        fprintf('Saving collected scenario predictions to %s', predTbPath);
        save(predTbPath, 'tbPred');
    end

    %% Sample weights

    logSampleWeight = NaN;
    if options.weights
        % compute importance weights based on baseline scenario with no reductions
        isBaseScenario = tbPredFull.("Catchment area load")=="1" &...
            tbPred.("Point load")=="1" & tbPred.("Internal load")=="1" & ...
            tbPred.("Atmospheric load") =="1";
        tbBaseScen = tbPredFull(isBaseScenario,:);
        if ~isempty(tbBaseScen)
            baseScenarioID = unique(tbBaseScen.scenarioID);
            if numel(baseScenarioID)>1
                warning('Multiple (%d) identifiers for base scenario',...
                    numel(baseScenario));
                baseScenarioID = baseScenarioID(1);
            end
            baseScenDir = fullfile(resultDir, ...
                sprintf('scenario_%s', baseScenarioID));
            sampleFile = fullfile(baseScenDir, 'sampleTable.mat');
            paramsFile = fullfile(baseScenDir, 'calibrationParams.mat');
            if ~ (isfile(sampleFile) && isfile(paramsFile))
                error(['Need both sample and parameter files, expecting:\n' ...
                    '%s\n%s'], sampleFile, paramsFile);
            end
            load(paramsFile, 'params');
            load(sampleFile, 'tbSamples');
            logSampleWeight = postPredImpWeights(tbBaseScen, tbSamples, params);
        else
            % use uniform weights
            warning('No predictions with base nutrient loads')
            logSampleWeight = 0;
        end
    else
        % uniform weights
        logSampleWeight = 0;
    end

    %% Summarising predictions

    tbChlaSum = scenario_summaries(tbPred, logSampleWeight, getGesTable(),...
        'quantiles', true);

    if options.saveSumTable
        sumTbPath = getNewSaveFile(...
            fullfile(resultDir, 'scenarioSummaryTable.mat'));
        save(sumTbPath, 'tbChlaSum');
        pltCmd = sprintf('plotChlaDistSets(''%s'')', sumTbPath);
        fprintf(['Saving scenario summary to %s.\n' ...
            'Summary can be plotted with %s.\n'], sumTbPath, pltCmd);
    end
    
    %% Plotting summary

    if options.plot
        plotChlaDistSets(tbChlaSum, "export",options.savePlot);
    end

end
