function [TT] = readOutput(outputpath, filtBlocks, options)
%READOUTPUT Reads FICOS simulator output into a timetable.
% TT = readOutput(outputDir) Reads FICOS simulator output from outputDir into a
% timetable.
%
% Optional arguments
% filtBlocks : (positional) list of FICOS's internal codes for water formations
% to load from the results. By default, all are loaded.
% 
% 'mode'      : In which mode to read the results, either 'ascii' or 'hdf5'
%              (default). 'ascii' mode is not recommended any more.
% 'vars'      : Load predictions for these variables.
% 'filtBlock' : The optional, positional filtBlocks argument can also be given
%               as a (name, value) argument.
% 
% Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    % Folder containing simulation results
    % either as a results.hdf5 file
    outputpath {mustBeFolder}
    % only load these blocks (simulation locations) from the results
    % if empty (default), load all blocks
    filtBlocks {mustBeA(filtBlocks, 'categorical')} = categorical()
    % 'hdf5'  : load results from results.hdf5 (default) or 
    % 'ascii' : from individual ascii files (deprecated)
    options.mode {mustBeTextScalar, ...
        mustBeMember(options.mode,{'ascii','hdf5'})} = 'hdf5';
    % which variables to load
    % load all by default
    options.vars {mustBeText} = getSimOutNames();
    % alternative way to define filtBlocks above
    options.filtBlocks {mustBeA(options.filtBlocks, 'categorical')} = ...
        categorical();
end

%% Ensuring given path is in the correct format
% Cast outputpath as char in case it was a string
outputpath = char(outputpath);
% Add '/' if missing
if(outputpath(end) ~= '/')
    outputpath = [outputpath '/'];
end
if(~isfolder(outputpath))
    error("Error: Folder %s does not exist.", outputpath)
end
%%

if isempty(filtBlocks) && ~isempty(options.filtBlocks)
    mustBeA(options.filtBlocks, 'categorical');
    filtBlocks = options.filtBlocks;
end

if strcmp(options.mode, 'ascii')
    % Load results from separate ascii files (legacy option)
    tcount = 0;
    timeout_max = 300;
    success = 0;
    
    varNameChars = alphanumericsPattern(1) | "_";
    varName = asManyOfPattern(varNameChars);
    asciiFilePattern = varName + "out.txt";
    
    TT = timetable();
    while(~success && tcount < timeout_max)
        try
            outFiles = dir(outputpath);
            % Select all ascii output files
            outFiles = outFiles(matches({outFiles(:).name}, asciiFilePattern));
            % Extract variable names from file names
            varNames = extractBefore({outFiles(:).name}, "out.txt");
            
            simData = struct();
            % Use first output as the template for output sizes 
            % and location columns
            [simData.(varNames{1}), polyids, mt] = loadAsciiPrediction(...
                outputpath, varNames{1});
            
            nTimes = numel(mt);
            % Measurement times
            % Simulation starts from t1 = 1/1/2006
            Aika = repmat(mt, numel(polyids),1);

            nVars = numel(varNames);
            % locations 
            polyID = repelem(int32(polyids(:)), numel(mt),1);

            simData.polyID = polyID;
            for i = 2:nVars
                simData.(varNames{i}) = ...
                    loadAsciiPrediction(outputpath, varNames{i});
            end
    
            TT = struct2table(simData);
            % Rename columns:
            newNames = cellfun(@(s)convertVarName(s), varNames, ...
                'UniformOutput', false);
            TT = renamevars(TT, varNames, newNames);
    
            TT = table2timetable(TT, 'RowTimes', Aika);
            TT.Properties.DimensionNames{1} = 'Aika';
            
            success = 1;
        catch
            pause(0.1);
            tcount = tcount + 0.1;
        end
        
    end
    
    if ~isempty(filtBlocks)
        TT = filterBlocks(TT, filtBlocks);
    end


elseif strcmp(options.mode, 'hdf5')
    % Read output directly from results.hdf5 (default)
    fpath = fullfile(outputpath, 'results.hdf5');
    % First load the time vector
    timeVec = h5read(fpath, '/Dates', [1 1], [3 Inf])';
    timeVec = datetime(timeVec);
    nt = numel(timeVec);


    % Get list of blocks
    blockInfo = h5info(fpath, '/Blocks');
    blockPaths = {blockInfo.Groups(:).Name}';

    if ~isempty(filtBlocks)
        % Filter to selected blocks
        blockPaths(~endsWith(blockPaths, string(filtBlocks))) = [];
    end
    wf_ids = categorical(extractAfter(blockPaths, '/Blocks/'));

    % Preallocate time table using first block as a template
    blkDatasets = h5info(fpath, blockPaths{1});
    blkDatasets = {blkDatasets.Datasets(:).Name}';
    if ~isempty(options.vars)
        blkDatasets(~matches(blkDatasets, options.vars)) = [];
    end
    TT = simOutToStruct(fpath, blockPaths{1}, blkDatasets, ...
        numel(timeVec));
  
    TT = struct2table(TT);
    TT = table2timetable(TT, "RowTimes", timeVec);
    
    nrows = nt*size(blockPaths,1);
    TT(nrows+1,:) = TT(1,:);
    TT(nrows+1,:) = [];

    for i = 2:size(blockPaths,1) 
        indT = (i-1)*nt + (1:nt);
        TT.block(indT) = wf_ids(i);
        TT.Time(indT) = timeVec;
        for j = 1:size(blkDatasets,1)
            TT.(blkDatasets{j})(indT) = ...
                h5read(fpath, ...
                    sprintf('%s/%s', blockPaths{i}, blkDatasets{j}),...
                    1, nt);
        end
    end

    if all(ismember(["cA"; "cC"], blkDatasets))
        % Estimate Chla biomass similarly to ascii results
        TT.chla = convertAFCtoChla(TT.Time, TT.cA, TT.cC);
    end
    
    % convert variable names to the ones used in likelihood function and
    % elsewhere
    oldVarNames = TT.Properties.VariableNames;
    newNames = cellfun(@(s)convertVarName(s), oldVarNames, ...
        'UniformOutput', false);
    TT = renamevars(TT, oldVarNames, newNames);
end

end
%% ------------------------------------------
function [predVec, polyIDs, timeVec] = loadAsciiPrediction(outputPath, varName)
    % Helper function for loading prediction from ASCII files
    fpath = fullfile(outputPath, [varName 'out.txt']); 
    fileCont = load(fpath);
    % Return as a long vector; first row is a header for locations
    predVec = fileCont;
    predVec(1,:) = [];
    predVec = predVec(:);

    % (Optional) return locations
    polyIDs = fileCont(1,:)';

    % (Optional) return time vector
    % first observations is on 1/1/2006
    timeVec = (datetime(2006,1,1) + ...
            caldays(0: (size(fileCont,1)-2)))';
end
% -------------------------
function newName = convertVarName(oldName)
    % Convert variable name to the one used in analysis

    switch(oldName)
        case {'ca_1_', 'cA'}
            newName = 'A';
        case {'cc_1_', 'cC'}
            newName = 'FC';
        case {'cdin_1_', 'cDIN_0'}
            newName = 'DIN1';
        case {'cdin_2_', 'cDIN_1'}
            newName = 'DIN2';
        case {'cdip_1_', 'cDIP_0'}
            newName = 'DIP1';
        case {'cdip_2_', 'cDIP_1'}
            newName = 'DIP2';
        case 'chla_1_'
            newName = 'chla';
        case 'rN2fixFC_'
            newName = 'rN2fixFC';
        otherwise
            % Remaining variables kept with the original name
            % Currently unused variables will have '_' at the end
            newName = oldName;
    end

end