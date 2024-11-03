function [tbPostPredQuant] = postPredIntervalBootstrap(F, y, w, ...
    obsModel, optSampling, optInterval)
    %POSTPREDINTERVALBOOTSTRAP Estimate posterior predictive credible 
    % intervals.
    %   tbPostPred = postPredIntervalBootstrap(F, y, w) Estimates posterior
    %   predictive credible intervals for posterior predictive simlations.
    %   
    %   Arguments:
    % 
    %   F  - (required) timetable containing posterior predictive simulations.
    %
    %   positional arguments:
    %
    %   y - observations used when sampling standard deviation parameter sigma 
    %   w - column vector of sample weights, one for each unique sampleID in F. 
    %
    %   (Name, value) arguments:
    %   - 'nBootstrap' - number of bootstrap samples
    %   used.
    %   - 'pInt' - credible interval width for simulation uncertainties
    %   - 'pObsInt' - credible interval width for posterior predictive 
    %     observation uncertainties
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        F (:,:) timetable
        y (:,:) timetable = load_intens(true, 'old', 'summer-algae');
        w (:,1) {mustBeNumeric} = 1;
        % structure defining observation model, i.e.
        % which variables are used, lower/upper bounds, transformations
        obsModel (1,1) struct = getObsModel();
        % how many bootstrap samples for the predictions
        optSampling.nBootstrap (1,1) {mustBePositive} = numel(w);
        % burn in used in Gibbs sampling for sigma
        optSampling.burnin (1,1) {mustBeNonnegative} = 100;
        % maximum amount of samples generated at once
        % 1 double = 8 bytes => 1e8 doubles ~= 800 Mb
        optSampling.maxChunk (1,1) {mustBePositive} = 1e7;
        % extra samples for sigmas or predictive observations
        optSampling.nSigmaSamp (1,1) {mustBePositive} = 1;
        optSampling.nYsamp (1,1) {mustBePositive} = 1;
        optInterval.intervalType (1,:) {mustBeTextScalar} = 'CI';
        optInterval.pInt (1,1) {mustBePositive, ...
            mustBeLessThanOrEqual(optInterval.pInt, 1)} = 0.95;
        optInterval.pObsInt (1,1) {mustBePositive, ...
            mustBeLessThanOrEqual(optInterval.pObsInt, 1)} = 0.80;
    end

    % Ensure simulation outputs and calibration data are in the same format
    y = cleanSimOutput(y);
    F = cleanSimOutput(F);

    %% Sigma samples
    fprintf('Sampling posterior sigmas.\n')
    tic
    [sigmaSamples, sigmaVariableGroups, sigmaSampleGroups] = ...
        samplePostPredSigmas( ...
        F, y, w, obsModel, ...
        nSamp=optSampling.nSigmaSamp, burnin=optSampling.burnin);
    toc
    clear('y')

    %% Extracting sample weights from tables if weights not given directly

    if isscalar(w)
        if ismember('logW', F.Properties.VariableNames)
            % take sample weights from F
            tbSamples = F(...
                F.wf_id == "1900000089" &...
                F.Time == datetime(2006,1,1) &...
                F.scenarioID=="1", {'sampleID','logW'});
            tbSamples = sortrows(tbSamples,'sampleID','ascend');
            w = exp(tbSamples.logW - logSumExp(tbSamples.logW));
        else
            % Use uniform weights if sample weights are not given
            w = ones(length(unique(F.sampleID)));
        end
    end

    %% Bootstrap sample for simulator parameterisations

    iSimSamp = datasample(1:numel(w), optSampling.nBootstrap, 'Weights',w);

    %% Preallocating result table

    F = stackVars(F, 'colName', 'prediction', 'idxName', 'variable');
    % Preallocate using predictions as template
    tbPostPredQuant = F;
    tbPostPredQuant(tbPostPredQuant.sampleID ~= categorical(1), :) = [];
    tbPostPredQuant = tbPostPredQuant(:, {'wf_id', 'variable', 'prediction'});

    quantCols = {'predLow', 'predMed', 'predHigh', 'obsLow','obsMed','obsHigh'};
    tbPostPredQuant{:, quantCols} = NaN;

    %%
    % When generating large amounts of samples, generate them in chunks to
    % avoid using too much memory at once
    totalRowSamp = optSampling.nBootstrap * optSampling.nSigmaSamp * ...
        optSampling.nYsamp;
    maxRows = floor(optSampling.maxChunk/totalRowSamp);
    % Using row size for 'DIN1' predictions for one sample as the template
    nRows = sum(F.variable=='DIN1' & F.sampleID=='1');
    nChunks = max(ceil(nRows/maxRows),1);


    %% Transform predictions before sampling
    % Observation model given for transformed observations
    F.prediction = max(F.prediction, 1e-2);
    F.prediction = obsModel.transFcn(F.prediction);


    %% Estimating intervals

    varNames = unique(sigmaVariableGroups);

    AFCseparately = all(ismember({'A','FC'}, sigmaVariableGroups));
    if AFCseparately
        varNames = setdiff(varNames, {'A','FC'});
    end


    cNames = categories(varNames);
    predVar = cell2struct(cNames, cNames);
    predSamp = structfun(@(v) extractPredSamples(F, v), predVar, ...
        'UniformOutput',false);
    % Corresponding time vector
    [~, tPred] = extractPredSamples(F, categorical("A"));

    predVar2 = predVar;
    if AFCseparately
        % Include chla in the resulting structures
        predVar2.chla = categorical("chla");
    end
    fInt = arrayfun( ...
        @(x) structfun(@(f) NaN(maxRows,3), predVar2, 'UniformOutput', false),...
        ones(nChunks,1));
    yInt = fInt;
    sSamp = structfun(@(v) extractSigmaSamples(...
        sigmaSamples, sigmaVariableGroups, v), predVar, ...
        'UniformOutput',false);
    for j = 1:nChunks
        fprintf(' chunk: %d/%d\n', j, nChunks);
        iChunkStart = (j-1)*maxRows+1;
        iChunkEnd = min(j*maxRows, nRows);
        chunkRows = iChunkStart:iChunkEnd;
        tChunk = tPred(chunkRows);

        predChunk = structfun(@(v) v(chunkRows,:), predSamp, ...
            'UniformOutput',false);
        % Intervals for simulations
        ySamp = structfun(@(v) estKernelSamples(v, iSimSamp, w, ...
            optSampling.nYsamp*optSampling.nSigmaSamp), ...
            predChunk, 'UniformOutput', false );
        ySamp = structfun(@(v) obsModel.retransFcn(v), ySamp, ...
            'UniformOutput',false);
        if(AFCseparately)
            % Chlorophyll a based on A and FC
            ySamp.chla = convertAFCtoChla(tChunk, ySamp.A, ySamp.FC);
        end
        fInt(j) = structfun(...
            @(y)estPredInterval(y, optInterval.intervalType, optInterval.pInt), ...
            ySamp, 'UniformOutput',false);

        fNames = cell2struct( fieldnames(predChunk), fieldnames(predChunk));
        % Intervals for observations
        ySamp = structfun(@(v) ...
            estPredSample(predChunk.(v), sSamp.(v), iSimSamp, optSampling.nYsamp), ...
            fNames, 'UniformOutput',false);
        ySamp = structfun(@(v) obsModel.retransFcn(v), ySamp, ...
            'UniformOutput',false);
        if AFCseparately
            ySamp.chla = convertAFCtoChla(tChunk, ySamp.A, ySamp.FC);
        end
        yInt(j) = structfun(...
            @(y)estPredInterval(y, optInterval.intervalType, ...
            optInterval.pObsInt), ySamp, 'UniformOutput',false);

    end
    % Combine intervals from all chunks into the timetable
    predVars = fieldnames(fInt);
    fCols = {'predLow', 'predMed', 'predHigh'};
    yCols = {'obsLow','obsMed','obsHigh'};
    for i = 1:numel(predVars)
        predVar = predVars{i};
        tbPostPredQuant{tbPostPredQuant.variable==predVar, fCols} = ...
            vertcat(fInt.(predVar));
        tbPostPredQuant{tbPostPredQuant.variable==predVar, yCols} = ...
            vertcat(yInt.(predVar));
    end

end



%% helper functions

function [Fsamp, tSamp] = extractPredSamples(tb, predVar)
    % Extract predictions for observations predVar for each parameter 
    % sample theta_i as a N x M matrix
    % where N is the length of a single prediction vector 
    % and M the number of parameter samples

    Fsamp = tb(tb.variable==predVar, {'wf_id', 'sampleID', 'prediction'});
    Fsamp = unstack(Fsamp, 'prediction', 'sampleID', ... 
        'GroupingVariables', {'Time','wf_id'});
    tSamp = Fsamp.Time;
    Fsamp = Fsamp{:, vartype("numeric")};
end

function Ssamp = extractSigmaSamples(sigmaSamples, gVar, predVar)
    % Extract sigma sample for a given predicted variable predVar
    % Output will be a n x m array where n is the number of sigma samples
    % m the number of simulator samples
    Ssamp = cell2mat( sigmaSamples(gVar==predVar)');
end

function ySamp = estPredSample(Fsamp, sigmaSamples, iFsamp, nSamp)
    % Generates posterior predictive samples given samples for mu and sigma
    Fsamp = repmat(Fsamp(:, iFsamp), 1, nSamp*size(sigmaSamples,1));
    sigmaSamples = sigmaSamples(:, iFsamp);
    sigmaSamples = repmat(sigmaSamples(:)', size(Fsamp,1), nSamp);
    ySamp = Fsamp + sigmaSamples .* normrnd(0,1, size(Fsamp));
end

function predInt = estPredInterval(ySamp, intType, pInterval)
    % Estimates posterior predictive interval for given samples
    switch intType
        case 'CI'
            qInterval = [(1-pInterval)/2, 0.5, (1+pInterval)/2];
            predInt = quantile(ySamp, qInterval, 2);
        otherwise
            error('Not implemented.')
    end

end

function ySamp = estKernelSamples(Fsamp, iFsamp, w, nSamp)
    % Kernel smoothed samples for simulator output

    % Using half of std for predictions, similarly to kernel1.m
    sigmaKernel = 0.5*std(Fsamp, w, 2);
    Fsamp = repmat(Fsamp(:, iFsamp), 1, nSamp);
    ySamp = Fsamp + diag(sigmaKernel) * normrnd(0,1, size(Fsamp));
end



