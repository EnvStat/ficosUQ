function [tbComp] = comparisonTable(F,y, obsModel, options)
%COMPARISONTABLE Join predictions and observations in a table.
%   tbComp = comparisonTable(F,y) Join predictions F with
%   observations y.
%
%   tbComp = comparisonTable(F,y,obsModel) Use the observation model obsModel
%   for the comparisons. If observation model contains transformations, the
%   resulting table will use transformed observations and predictions.
%
%   Additional named arguments:
%   'predVar' : (text scalar) Use this name to as the column name for predictions 
%               (default 'prediction').
%   'obsVar'  : (text scalar) Use column name for observations (default
%               'observation'.
%   'markCensored' : adds a column denoting whether an observation is below the
%                    censoring lower bound.
%   'truncateToLb' : Truncate values to lower bounds. (default 'off')
%   'idxName'      : name of the indexing variable when pivoting (stacking) the
%                    table into the tall format (default: 'variable').
%
%   See also getObsModel.
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    F (:,:) timetable
    y (:,:) timetable
    obsModel (1,1) struct = getObsModel();
    options.truncateToLb {mustBeTextScalar, mustBeMember(...
        options.truncateToLb, {'off', 'on', 'simLB', 'censureLB'})} = 'off';
    options.markCensored {mustBeNumericOrLogical} = false;
    options.predVar {mustBeTextScalar} = 'prediction';
    options.obsVar {mustBeTextScalar} = 'observation';
    options.idxName {mustBeTextScalar} = 'variable';
end

F = cleanSimOutput(F);
y = cleanSimOutput(y);

% Pivot into stacked format for observed variables
idxName = options.idxName;
predVar = options.predVar;
obsVar = options.obsVar;

if ~ismember(idxName, F.Properties.VariableNames)
    F = stackVars(F, obsModel.obsVars, 'colName', predVar, ...
        'idxName', idxName);
end
if ~ismember(idxName, y.Properties.VariableNames)
    y = stackVars(y, obsModel.obsVars, 'colName', obsVar, ...
        'idxName', idxName);
end

% Select only relevant columns
sharedCols = {'wf_id', idxName};
fCols = [sharedCols, {predVar, 'sampleID'}];
fCols = intersect(fCols, F.Properties.VariableNames);
F = F(:, fCols);
y = y(:, [sharedCols, {obsVar}]);

% Match observations with predictions
tbComp = innerjoin(y, F, 'Keys', {'Time', 'wf_id', idxName});


% Transform based on observation model
tbComp.(predVar) = obsModel.transFcn(tbComp.(predVar));
tbComp.(obsVar) = obsModel.transFcn(tbComp.(obsVar));

% Add columns for transformed lower bounds 
tbLB = table(categorical(obsModel.obsVars(:)), ...
    obsModel.transFcn(obsModel.censureLB(:)), ...
    obsModel.transFcn(obsModel.simLB(:)), ...
    'VariableNames',{idxName, 'censureLB', 'simLB'});
tbComp = innerjoin(tbComp, tbLB, 'Keys', {idxName});

if options.markCensored
    tbComp.censored = tbComp.(obsVar) < tbComp.censureLB;
end

switch options.truncateToLb
    case 'simLB'
        tbComp.(predVar) = max(tbComp.(predVar), tbComp.simLB);
        tbComp.(obsVar) = max(tbComp.(obsVar), tbComp.simLB);
    case 'censureLB'
        tbComp.(predVar) = max(tbComp.(predVar), tbComp.censureLB);
        tbComp.(obsVar) = max(tbComp.(obsVar), tbComp.censureLB);
    case 'on'
        tbComp.(predVar) = max(tbComp.(predVar), tbComp.simLB);
        tbComp.(obsVar) = max(tbComp.(obsVar), tbComp.censureLB);
end


end

