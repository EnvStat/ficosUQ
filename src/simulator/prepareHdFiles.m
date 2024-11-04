function [targetFilePath] = prepareHdFiles(options)
    %PREPAREHDFILES Prepares the hd_files.hdf5 input file used by the 
    % FICOS simulator.
    %   hdFilePath = prepareHdFiles(...) Internal function used for preparing
    %   the boundary condition file hd_files.hdf5 for the FICOS simulator.
    % 
    %   Wrapper for the 'solar_ma.py' Python script.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.hdf5dir {mustBeTextScalar} = 'data/';
        options.baseFileName {mustBeTextScalar} = 'hd_files.hdf5';
        options.targetFileName {mustBeTextScalar} = '';
        options.movingAvgFileNameFormat {mustBeTextScalar} = ...
            'hd_files_ma_%d.hdf5';
        options.numDays (1,1) {mustBeNonnegative} = 0;
        options.pythonDir {mustBeFolder} = fullfile('src','python/');
        % python3 interpeter command
        % use virtual environment by default
        options.pythonInterpreter {mustBeTextScalar} = ...
            getPythonInterpreter('fullPath',true);
        options.pythonScriptName {mustBeTextScalar} = 'solar_ma.py';
        options.targetDir {mustBeFolder} = 'data/';
    end

    baseFilePath = fullfile(options.hdf5dir, options.baseFileName);
    if isempty(options.targetFileName)
        if options.numDays
            options.targetFileName = ...
                sprintf(options.movingAvgFileNameFormat, options.numDays);
        else
            options.targetFileName = options.baseFileName;
        end
    end
    if options.numDays == 0
        % no need to compute moving average, use the base file
        sourceFilePath = fullfile(options.hdf5dir, options.baseFileName);
    else
        % See if there is already a file with moving average computed previously
        sourceFileName = sprintf(options.movingAvgFileNameFormat,...
            options.numDays);
        sourceFilePath = fullfile(options.hdf5dir, sourceFileName);
    
        
        
        % store current working directory
        pwdOrig = pwd;
        if ~isfile(sourceFilePath)
            % create file first with a python script
    
            % script assumes input file is in the current directory
            % operating in script directory, then moving the result
            copyfile(baseFilePath, options.pythonDir);
            cd(options.pythonDir);
            tmpBaseFile = fullfile('.', options.baseFileName);
            tmpSrcFile = fullfile(options.pythonDir, sourceFileName);
            % 
            pythonScriptPath = options.pythonScriptName;
            pythonCmd = sprintf('%s %s %d', ...
                options.pythonInterpreter, pythonScriptPath, options.numDays);
            [status, cmdOut] = system(pythonCmd);
            % delete temporary source file and return to original folder
            delete(tmpBaseFile);
            cd(pwdOrig);
        
            if status
                if isfile(tmpSrcFile)
                    delete(tmpSrcFile)
                end
                error(cmdOut);
            end
            % move result to hdf5 folder, so it can be re-used next time
            movefile(tmpSrcFile, sourceFilePath);       
        end
    end
    % copy target file to target folder
    targetFilePath = fullfile(options.targetDir, options.targetFileName);
    if strcmp(targetFilePath, baseFilePath)
        [~, tgtFileName, fExt] = fileparts(options.targetFileName);
        newTargetFileName = sprintf('%s_%s%s', ...
            tgtFileName, getTimeStamp(2), fExt);
        targetFilePath = fullfile(options.targetDir, newTargetFileName);

        warning(['Attempting to overwrite base hdf5-file: %s. ' ...
            'Writing instead to %s'], baseFilePath, targetFilePath);
    end

    if strcmp(targetFilePath, sourceFilePath)
        % file already exists
        return
    end
    if isfile(targetFilePath)
        % delete previous file to ensure copy is used
        delete(targetFilePath);
    end
    copyfile(sourceFilePath, targetFilePath);
end

