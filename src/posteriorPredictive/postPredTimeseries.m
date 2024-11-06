function [tbPostPredSum, tbSamples] = postPredTimeseries(resultDir, options)
    %POSTPREDTIMESERIES Posterior predictive time series
    %   tbPostPred = postPredTimeseries(resultDir) Summarise posterior
    %   predictive timeseries from simulations within resultDir.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    arguments
        resultDir {mustBeFolder}
        options.useImportanceWeights {mustBeNumericOrLogical} = true;
        options.observations {mustBeA(options.observations, 'timetable')} = ...
            load_intens(true, 'old', 'summer-algae');
        options.obsModel struct = getObsModel();
        options.computePostDens {mustBeNumericOrLogical} = true;
        options.plot {mustBeNumericOrLogical} = true;
        options.savePlot {mustBeNumericOrLogical} = false;
        options.saveIntermResults {mustBeNumericOrLogical} = true;
        options.showWarnings (1,1) {mustBeNumericOrLogical} = false;
    end

    if ~options.showWarnings
        warnStruct = warning;
        warning('off', 'MATLAB:table:ModifiedVarnamesUnstack');
    end

    %% Collect raw results
    

    tbPostPred = collectPostPredResults(resultDir, ...
        'vars', getSimOutNames('calibration'), ...
        'filtBlocks',filterBlocksArg('intens'));
    
    % Load additional information about parameters if able
    tbSamples = table();
    sampleFile = fullfile(resultDir, 'sampleTable.mat');
    paramsFile = fullfile(resultDir, 'calibrationParams.mat');
    if isfile(sampleFile)
        load(sampleFile, 'tbSamples');
        samplesLoaded = true;
    else
        samplesLoaded = false;
    end
    
    if isfile(paramsFile)
        load(paramsFile, 'params');
        paramsLoaded = true;
    else
        paramsLoaded = false;
    end

    %% Sample weight information  
    uqSamp = unique(tbPostPred.sampleID);
    nUqSamp = length(uqSamp);
    if options.useImportanceWeights
        if samplesLoaded && paramsLoaded && options.computePostDens
            % compute importance weights
            tbSamples = tbSamples(ismember(tbSamples.sampleID, uqSamp),:);
            logP = NaN(nUqSamp,1);
            for i = 1:nUqSamp
                idSamp = uqSamp(i);
                xSamp = tbSamples.sample(tbSamples.sampleID==idSamp,:);
                fSamp = tbPostPred(tbPostPred.sampleID==idSamp,:);
                fSamp = cleanSimOutput(fSamp, 'old');
                logP(i) = -1*params.post_f(xSamp, fSamp);
            end
            tmp = table(uqSamp(:), logP, 'VariableNames',{'sampleID', 'logP'});
            tbSamples.logP = [];
            tbSamples.logW = [];
            tbSamples = innerjoin(tbSamples, tmp, 'Keys', {'sampleID'});
            tbSamples.logW = tbSamples.logP - tbSamples.logQ;
        elseif samplesLoaded && all(~isnan(tbSamples.logW))
            % use already defined weights
            tbSamples = tbSamples(ismember(tbSamples.sampleID, uqSamp),:);
        else
            % uniform weights
            tbSamples = table(uqSamp(:), zeros(nUqSamp,1),...
                'VariableNames',{'sampleID', 'logW'});
        end
        % Join sample information with predictions
        tmp = tbSamples(:,{'sampleID', 'logW'});
        tbPostPred = innerjoin(tbPostPred, tmp, 'Keys', {'sampleID'});
        % weights in linear scale
        w = tbSamples.logW(:);
        w = exp(w - logSumExp(w));
    else
        % use uniform weights
        w = ones(nUqSamp,1);
        if ismember('logW', tbPostPred.Properties.VariableNames)
            tbPostPred.logW = [];
        end
    end
    %% Posterior predictive summary

    tbPostPredSum = postPredIntervalBootstrap(tbPostPred, ...
        options.observations, w,...
        options.obsModel);

    if options.saveIntermResults
        resultPath = fullfile(resultDir, 'postPredSummary.mat');
        if isfile(resultPath)
            % save to new file with timestamp if file already exists
            [~, resFile, resExt] = fileparts(resultPath);
            resFileNew = sprintf('%s_%s%s', resFile, getTimeStamp(), resExt);
            resultPath = fullfile(resultDir, resFileNew);
        end
        fprintf('Saving summary table to %s\n', resultPath);
        pltFun = 'plotPostPredTimeseries';
        fprintf(['Results can be plotted with: \n' ...
            '%s(''%s'')'], pltFun, resultPath);
        save(resultPath, 'tbPostPredSum');
    end

    %% Plotting
    if options.plot
        plotPostPredTimeseries(tbPostPredSum, 'save', options.savePlot);
    end


    if ~options.showWarnings
        warning(warnStruct);
    end
end