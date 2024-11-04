function [pyInterp] = getPythonInterpreter(options)
    %GETPYTHONINTERPETER Get path to Python interpreter.
    %   pyInterp = getPythonInterpreter() Returns default path to python3
    %   interpreter needed for running this ficosUQ project's Python scripts.
    %
    %   Used by internal functions.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.relInterpPath = 'src/python/pyEnv/venv/bin/python3';
        % require full path, otherwise return relative path
        options.fullPath {mustBeNumericOrLogical} = false;
    end

    if options.fullPath
        [pwdRoot, pwdName, ~] = fileparts(pwd);
        pyInterp = fullfile(pwdRoot, pwdName, options.relInterpPath);
        if ~isfile(pyInterp)
            error(['getPythonInterpreter:Invalid path: %s.\n' ...
                'Make sure that python libaries are installed and ' ...
                'that you are using this package from the project root folder.'], ...
                pyInterp)
        end
    else
        pyInterp = options.relInterpPath;
    end
    
end

