function [depOut] = findDependencies(sourceFiles, options)
    %FINDDEPENDENCIES Lists dependencies for this projects Matlab source files.
    % Find dependencies for given matlab sourceFiles, i.e.
    % which Matlab Toolboxes and which Matlab sourcefiles are needed to
    % run them.
    %   Example: Find which source local source files are needed for the main
    % calibration script as follows.
    % >> string(getfield(findDependencies("mode",'local'),'files'))
    %
    % If the ficosUQ
    %
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        sourceFiles {mustBeFile} = {fullfile('src/ficos_calibration.m')};
        % structure for the output
        % 'raw' - simply return the default output without any other structure
        % '' - for files, return local paths
        options.mode {mustBeTextScalar, ...
            mustBeMember(options.mode, {'raw', 'local'})} = 'raw';
        options.projectRoot {mustBeTextScalar} = 'ficosUQ/';
        options.includeSubmodules {mustBeNumericOrLogical} = false;
        % divide file list output by directory:
        options.subStructStyle {mustBeTextScalar, mustBeMember(...
            options.subStructStyle, {'none','mainFolder'})} = ...
            'none';
        % include this function
        options.includeSelf {mustBeNumericOrLogical} = true;
    end

    currPath = pwd;
    [~, currDir] = fileparts(currPath);
    srcDir = fullfile(options.projectRoot, 'src');
    if ~isfolder(srcDir)
        error("Expecting to find folder %s. " + ...
            "Try running this function with the named argument" + ...
            "'projectRoot',''%s''.", currDir);
    end

    if options.includeSelf
        sourceFiles = union(sourceFiles, {'src/util/findDependencies.m'});
    end
    [fileList, prodList] = matlab.codetools.requiredFilesAndProducts(...
        sourceFiles);

    depOut = struct;
    switch options.mode
        case 'raw'
            depOut.files = fileList(:); % easier to read vertically
            depOut.products = prodList;
        case 'local'
            % Return Matlab products (toolboxes) as is
            depOut.product = prodList;
            % Additional structure for files:

            % begin file paths from project dir
            fileList = extractAfter(fileList(:), options.projectRoot);


            % handle files in project root separately
            inSourceDir = startsWith(fileList, 'src/');
            rootFiles = fileList(~inSourceDir);

            sourceFiles = extractAfter(fileList(inSourceDir), 'src/');
            % similarly handle submodules separately
            isSubModule = startsWith(sourceFiles, 'submodules/');
            subModFiles = sourceFiles(isSubModule);
            sourceFiles = sourceFiles(~isSubModule);

            if options.includeSubmodules
                subModFiles = extractAfter(subModFiles, 'submodules/');
            else
                % discard
                subModFiles = {};
            end

            switch options.subStructStyle
                case 'none'
                    % simply return the remaining paths
                    depOut.files = [rootFiles(:); sourceFiles(:); subModFiles(:)];
                case 'mainFolder'
                    % structures for each main folder
                    depOut.files.root = rootFiles(:);
                    %
                    depOut.files.source = sourceFiles(:);
                    %
                    if options.includeSubmodules
                        depOut.files.submodules = isSubModule;
                    end
            end
    end
end