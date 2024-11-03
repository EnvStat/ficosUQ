function [pyInterp] = getPythonInterpeter(options)
    %GETPYTHONINTERPETER Get path to Python interpreter.
    %   pyInterp = getPythonInterpreter() Returns default path to python3
    %   interpreter needed for running this ficosUQ project's Python scripts.
    %
    %   Used by internal functions.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.projectRootName = 'ficosUQ';
        options.relInterpPath = 'src/python/pyEnv/venv/bin/python3';
        % require full path, otherwise return relative path
        options.fullPath {mustBeNumericOrLogical} = false;
    end

    if options.fullPath
        [pwdRoot, pwdName, ~] = fileparts(pwd);
        if ~strcmp(pwdName, options.projectRootName)
            error(['getPythonInterpeter:Current working directory needs ' ...
                'to be at project root. Expecting: %s, path is instead: %s'], ...
                options.projectRootName, pwd);
        end
        pyInterp = fullfile(pwdRoot, pwdName, options.relInterpPath);
    else
        pyInterp = options.relInterpPath;
    end
    
end

