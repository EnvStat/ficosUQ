function [num_cores] = getNumCores(def_cores, options)
    %GETNUMCORES Retrieves the number of active or available cpu cores. 
    %   num_cores = getNumCores() Retrieves the number of cpu cores actively 
    %   in the current parallel pool. If there are no current parallel pools, 
    %   instead return the amount of available cores.
    %
    %   num_cores = getNumCores(...,'def_cores', n) Return n as a fall back
    %   value in case the number of cores could not be determined.
    %
    %   num_cores = getNumCores(...,'poolobj', poolobj) Determine number of 
    %   cores from poolobj instead of the current parallel pool.
    %
    %   num_cores = getNumCores(...,'showWarning',true) Show warnings if 
    %   retrieval fails.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments 
        def_cores (1,1) {mustBePositive, mustBeInteger} = 20;
        % ---- additional options -----
        options.poolobj = gcp('nocreate');
        % whether to show warning if retrieval fails 
        options.showWarning {mustBeNumericOrLogical} = false;
    end

    try
        if ~isempty(options.poolobj) && isa(options.poolobj, 'parallel.Pool')
            % Use the number of active workers
            num_cores = options.poolobj.NumWorkers;
        else
            % Use undocumented MATLAB feature to get number of available cores
            % https://stackoverflow.com/questions/8311426/how-can-i-query-the-number-of-physical-cores-from-matlab
            num_cores = feature('numcores');
        end
    catch ME
        if options.showWarning
            warning(ME);
        end
        num_cores = def_cores;
    end
end

