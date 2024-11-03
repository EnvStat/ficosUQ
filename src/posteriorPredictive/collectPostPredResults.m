function [postPredResult] = collectPostPredResults(resultRootDir, ...
        outputOpt, options, optMtdt)
    %COLLECTPOSTPREDRESULTS Collect posterior predictive simulation results.
    %   postPredResults = collectPostPredResults(resultRootDir) Collect
    %   posterior predictive simulation results from within resultRootDir.
    %
    %   ... = collectPostPredResults(dir, 'recurse', true) Collect recursively
    %   from subfolders within dir. Used for nutrient load scenario predictions.
    %   Only recurses only level.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        resultRootDir {mustBeFolder}
        % options to readOutput
        outputOpt.mode {mustBeTextScalar} = 'hdf5';
        outputOpt.vars = getSimOutNames('extScenario');
        outputOpt.filtBlocks ...
            {mustBeA(outputOpt.filtBlocks, 'categorical')} = categorical();
        % other options
        options.recurse {mustBeNumericOrLogical} = false;
        % post process results
        options.postFcn {mustBeScalarOrEmpty} = [];
        % metadata formatting options
        % names used in encoded metadata
        optMtdt.mtdtNames = ...
            {'catchmentLoad', 'pointLoad', 'internalLoad', 'atmosphericLoad'};
        % names used for summary tables
        optMtdt.mtdtNamesTb = ...
            {'Catchment area load','Point load','Internal load', 'Atmospheric load'};
    end
    
    subdirs = dir(resultRootDir);
    subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));
    subdirs = subdirs([subdirs.isdir]);
    ndirs = length(subdirs);

    % use first results as template
    if options.recurse
        % go into subsub dirs
        % used for scenario predictions
        subPathFcn = @(s) fullfile(s.folder, s.name);
        outputArgs = [namedargs2cell(outputOpt), namedargs2cell(optMtdt)];
        collectFcn = @(x) collectPostPredResults(subPathFcn(x), ...
            outputArgs{:}, 'recurse', false);
    else
        collectFcn = @(x) collectResult(x, outputOpt, optMtdt);
    end
    
    if ~isempty(options.postFcn) && isa(options.postFcn, 'function_handle')
        resFcn = @(x)options.postFcn(collectFcn(x));
    else
        resFcn = collectFcn;
    end

    postPredResult = resFcn(subdirs(1));
    nPred = size(postPredResult,1);
    timeVar = postPredResult.Properties.DimensionNames{1};
    tVec = postPredResult.(timeVar);
    nTotal = ndirs*nPred;
    postPredResult(nTotal+1,:)= postPredResult(1,:);
    postPredResult(nTotal+1,:)= [];
    
    for i = 2:ndirs
        istart = (i-1)*nPred+1;
        iend = i*nPred;
        ind = (istart:iend);
        postPredResult(ind,:) =  resFcn(subdirs(i));
        postPredResult.(timeVar)(ind) = tVec;
    end

end

%% ------------

function [tbPred] = collectResult(subdir, outputOpt, optMtdt)
    % fetch one result
    resultDir = fullfile(subdir.folder, subdir.name);
    resultFile = fullfile(resultDir, 'results.hdf5');
    
    % Metadata
    resMtdt = readScenarioMetadata(resultFile);
    mtdtVars = fieldnames(resMtdt);
    resMtdt = struct2table(resMtdt);
    resMtdt = convertvars(resMtdt, mtdtVars, 'categorical');
    resMtdt = renamevars(resMtdt, ...
        optMtdt.mtdtNames, optMtdt.mtdtNamesTb);

    % predictions
    outputArgs = namedargs2cell(outputOpt);
    tbPred = readOutput(resultDir, outputArgs{:});

    % add metadata
    nPred = size(tbPred, 1);
    resMtdt = repmat(resMtdt, nPred, 1);
    tbPred = [tbPred, resMtdt];
end