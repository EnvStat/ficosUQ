function [gpstuffroot] = getGpstuffroot(options)
    %GETGPSTUFFROOT Retrieve path to gpstuff installation root.
    %   gpstuffroot = getGpstuffroot(); Retrieves gpstuffroot relative to the
    %   project root folder.
    %
    %   gpstuffroot = getGpstuffroot('absolute', true) Returns absolute path
    %   to gpstuff installation. Has to be run from the project root folder.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.absolute (1,1) {mustBeNumericOrLogical} = false;
        % paths defined as arguments here in case they are changed later
        options.relGpstuffPath {mustBeTextScalar} = ...
            fullfile('src','submodules/','gpstuff/');
    end

    if ~options.absolute
        gpstuffroot = options.relGpstuffPath;
    else
        [parentDir, currDir] = fileparts(pwd);
        gpstuffroot = fullfile(parentDir, currDir, options.relGpstuffPath);
        if ~isfolder(gpstuffroot)
            error(['getGpstuffroot:Invalid path: %s.\n ' ...
                'Run this function from the project root folder, ' ...
                'currently in: %s'], gpstuffroot, fullfile(parentDir,currDir));
        end
    end
end

