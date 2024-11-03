function [gpstuffroot] = getGpstuffroot(options)
    %GETGPSTUFFROOT Retrieve path to gpstuff installation root.
    %   gpstuffroot = getGpstuffroot(); Retrieves gpstuffroot relative to the
    %   ficosUQ folder.
    %
    %   gpstuffroot = getGpstuffroot('absolute', true) Returns absolute path
    %   to gpstuff installation. Has to be run from the ficosUQ folder.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.absolute (1,1) {mustBeNumericOrLogical} = false;
        % paths defined as arguments here in case they are changed later
        options.relGpstuffPath {mustBeTextScalar} = ...
            fullfile('src','submodules/','gpstuff/');
        options.projectRootName = 'ficosUQ';
    end
    
    if ~options.absolute
        gpstuffroot = options.relGpstuffPath;
    else
        [parentDir, currDir] = fileparts(pwd);
        if ~strcmp(currDir, options.projectRootName)
            error('Must be in the %s project folder. Currently in: %s', ...
                options.projectRootName, pwd);
        end
        gpstuffroot = fullfile(parentDir, currDir, options.relGpstuffPath);
    end
end

