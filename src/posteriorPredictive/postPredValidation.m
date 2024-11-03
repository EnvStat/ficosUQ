function [comparisonStruct] = postPredValidation(options)
    % POSTPREDVALIDATION Compute model performance metrics for posterior
    % predictions.
    %
    % vMetrics = postPredValidation('resultFile', postPredSumFile) Calculates a
    % set of model performance metrics for posterior predictive summaries saved
    % in given file.
    %
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.resultFile = fullfile('results/calibration_example/postPredSummary.mat');
        options.calibrationData = 'default';
        options.transform ...
            {mustBeTextScalar, ...
            mustBeMember(options.transform, {'log', 'sqrt', 'none'})} = 'log';
        options.censoring (1,1) {mustBeNumericOrLogical} = true;
        options.moreMetrics (1,1) {mustBeNumericOrLogical} = true;
        % Save all metrics to one table (true) or to separate tables (false)
        options.CombinedTable (1,1) {mustBeNumericOrLogical} = true;
        options.digits (1,1) {mustBeNonnegative} = 2;
        options.Extension = 'xlsx';
    end
    
    %% load predictive medians
    load(options.resultFile, ...
        'tbPostPredSum');
    
    tbPred = tbPostPredSum(:, {'wf_id', 'variable', 'predMed'});

    
    %% Load calibration data
    
    if(strcmp(options.calibrationData, 'default'))
        tt_data = load_intens(true, 'old', 'summer');
        tbObsCalib = cleanSimOutput(tt_data);
    else
        error('Invalid calibration data')
    end
    
    %% Construct comparison table
    
    tbPred.Properties.DimensionNames{1} = 'Aika';
    tbObsCalib.Properties.DimensionNames{1} = 'Aika';
    
    truncateOpt = 'off';
    switch options.transform
        case 'log'
            transFcn = @(x) log(x);
            truncateOpt = 'simLB';
        case 'sqrt'
            transFcn = @(x) sqrt(x);
        case 'none'
            transFcn = @(x) x;
    end
    oM = getObsModel("transFcn",transFcn);
    if ~ismember('chla', oM.obsVars)
        oM.obsVars = [oM.obsVars, {'chla'}];
        oM.censureLB = [oM.censureLB, 0.01];
        oM.simLB = [oM.simLB, 0.01];
    end

    tbComp = comparisonTable(tbPostPredSum, tt_data, oM, ...
        'truncateToLb',truncateOpt, 'markCensored',options.censoring, ...
        'predVar','predMed');
    if options.censoring
        % only include non censored observation for performance metrics
        tbComp(tbComp.censored,:) = [];
    end
    
    %% Compute model performance metrics
       
    locations = unique(tbComp.wf_id);
    varNames = unique(tbComp.variable);
    
    tmp = NaN(numel(locations), numel(varNames));
    tbR2 = array2table( tmp, ...
        'RowNames',string(locations), ...
        'VariableNames', string(varNames));
    
    comparisonStruct = struct();
    if options.moreMetrics
        % Other criteria similarly
        comparisonStruct.R2 = tbR2;
        comparisonStruct.RMSE = tbR2;
        % RSR : RMSE/SD(obs)
        comparisonStruct.RSR = tbR2;
        % Corr(obs, prediction)
        comparisonStruct.corr = tbR2;
        % Correlation squared
        comparisonStruct.corr2 = tbR2;
        % Index of agreement
        comparisonStruct.IndAgr = tbR2;
        % Percent Bias : sum_i(obs_i - pred_i)/sum_i(obs_i) 
        comparisonStruct.PBIAS = tbR2;
        % Mean avaerage error : sum_i |obs_i - pred_i|/n
        comparisonStruct.MAE = tbR2;

    end
    
    for i = 1:length(locations)
        for j = 1:length(varNames)
            tbTmp = tbComp(tbComp.wf_id == locations(i) & ...
                tbComp.variable == varNames(j),:);
            if options.moreMetrics
                comparisonStruct.R2{i,j} = ...
                    r2Fcn(tbTmp.observation, tbTmp.predMed);
                comparisonStruct.RMSE{i,j} = ...
                    rmse(tbTmp.predMed, tbTmp.observation);
                comparisonStruct.RSR{i,j} = ...
                    rsrFcn(tbTmp.observation, tbTmp.predMed);
                comparisonStruct.corr{i,j} = ...
                    corr(tbTmp.observation, tbTmp.predMed);
                % correlation squared
                comparisonStruct.corr2{i,j} = comparisonStruct.corr{i,j}.^2;
                comparisonStruct.IndAgr{i,j} = ...
                    indAgreeFcn(tbTmp.observation, tbTmp.predMed);
                comparisonStruct.PBIAS{i,j} = ...
                    pbiasFcn(tbTmp.observation, tbTmp.predMed);
                comparisonStruct.MAE{i,j} = ...
                    maeFcn(tbTmp.observation, tbTmp.predMed);
            else
                tbR2{i,j} = r2Fcn(tbTmp.observation, tbTmp.predMed);
            end

        end
    end
    
    % convert location names
    tbR2.Properties.RowNames = wf_id_to_station(tbR2.Properties.RowNames);
    
    %% Export
    resultDir = fullfile('results/validation/', getTimeStamp());
    mkdir(resultDir);

    if options.moreMetrics
        if options.CombinedTable
            % comparison 
            tbRes = combineTables(comparisonStruct);
            % rounding
            tbRes = roundFcn(tbRes, options.digits);
            resultFileStem = 'metrics';

            if strcmp(options.transform, 'none')
                resultFileName = sprintf('%s_original_scale.%s', ...
                    resultFileStem, options.Extension);
            else
                resultFileName = sprintf('%s_%s_scale.%s', ...
                    resultFileStem, options.transform , options.Extension);
            end
            writetable(tbRes, fullfile(resultDir,resultFileName),...
                'WriteRowNames',false);
        else
            fieldNames = fieldnames(comparisonStruct);
            resDir = fullfile('results/validation', getTimeStamp());
            mkdir(resDir);
            for ii = 1:length(fieldNames)
                tbRes = comparisonStruct.(fieldNames{ii});
                % Round values
                tbRes.Variables = round(tbRes.Variables, options.digits);
                exportFcn(tbRes, ...
                    fieldNames{ii}, resDir);
            end
        end

    else
        resultFileName = 'R2.csv'; 
        writetable(tbR2, fullfile(resultDir, resultFileName), 'WriteRowNames',true);
    end

end
%% ---------------
function [r2] = r2Fcn(obs, pred)
    % Coefficient of determination
    SSR = sum((obs - pred).^2);
    SST = sum(obs.^2);
    r2 = 1 - SSR/SST;
end

function [dAgr] = indAgreeFcn(y, f)
    % Index of agreement
    % Using notation in (Moriasi et al. 2015)
    sqsum = sum( (y-f).^2);
    my = mean(y);
    devPred = abs(f -my);
    devObs = abs(y - my);
    devSqsum = sum( (devPred + devObs).^2);
    dAgr = 1 - sqsum/devSqsum;
end
% ----
function [RSR] = rsrFcn(y, f)
    % RSR: RMSE to observations Standard deviation Ratio
    nom = sqrt(sum( (y - f).^2));
    denom = sqrt(sum( (y - mean(f)).^2));
    RSR = nom/denom;
end
% ----
function [pB] = pbiasFcn(y, f)
    % Percent bias
    biasSum = sum(y - f);
    obsSum = sum(y);
    pB = 100*biasSum/obsSum;
end
% ----
function [MAE] = maeFcn(y, f)
    % Mean average error
    MAE = mean(abs(y - f));
end
% ----
function [outArg] = exportFcn(tbRes, metricName, resultDir)
    % export each comparison metric
    arguments
        tbRes {mustBeA(tbRes, {'table','timetable'})}
        metricName {mustBeTextScalar}
        resultDir {mustBeFolder} = 'results/validation';
    end
    % convert location names
    tbRes.Properties.RowNames = wf_id_to_station(tbRes.Properties.RowNames);
 
    if strcmp(resultDir, 'results/validation')
        % add timestamp to filename if timestamp is not in the folder
        resultFileName = sprintf('%s_%s.csv', getTimeStamp(), metricName);
    else
        % otherwise filename is just the name of the metric
        resultFileName = sprintf('%s.csv',metricName);
    end
    writetable(tbRes, fullfile(resultDir, resultFileName), ...
        'WriteRowNames',true);
    outArg = 0; % success
end
% ------
function [tbOut] = combineTables(compStruct, colName)
    arguments
        compStruct struct
        colName {mustBeTextScalar} = 'Observable';
    end
    % combine all result tables into one table

    % Transpose each table
    compStruct = structfun(@(tb) transposeTb(tb, colName), ...
        compStruct, UniformOutput=false);
    % convert to table
    tbOut = struct2table(compStruct);
    % split nested structure
    tbOut = splitvars(tbOut);
    
    % Trim duplicate columns
    colNames = tbOut.Properties.VariableNames;
    dupCols = colNames(endsWith(colNames, colName));
    % Use first of these for the actual variable
    tbOut = renamevars(tbOut, dupCols(1), {colName});
    % Discard rest
    tbOut = removevars(tbOut, dupCols(2:end));
end
% -----
function [tbOut] = transposeTb(tb, colName)
    % helper function for tranposing result tables
    arguments
        tb {mustBeA(tb, 'table')}
        colName {mustBeTextScalar} = 'observable';
    end
    % helper variable for stations
    tb.Station = wf_id_to_station(tb.Properties.RowNames);

    tbOut = rows2vars(tb, "VariableNamingRule","preserve", ...
        "VariableNamesSource", "Station");
    
    tbOut.OriginalVariableNames = categorical(tbOut.OriginalVariableNames);
    tbOut = renamevars(tbOut, 'OriginalVariableNames', colName);
end
% ----------
function [tb] = roundFcn(tb, digits)
    % helper function for rounding all numeric table columns
    tb(:, vartype("numeric")) = ...
        round( tb(:, vartype("numeric")), digits);
end