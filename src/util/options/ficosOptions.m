function [ficosOpt] = ficosOptions(ficosOpt, optFields)
    %FICOSOPTIONS Set options for running FICOS simulator.
    %   ficosOpt = ficosOptions() Constructs a structure with default paths 
    %   for running FICOS.
    %
    %   ficosOpt = ficosOptions('hdf5Opt', hdf5options()) Include also default
    %   hdf5 paths.
    % 
    %   See also RUN_SIM
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        ficosOpt.simPath {mustBeTextScalar} = fullfile('bin/ficos');
        % file names for simulator and its .ini file
        ficosOpt.simFileName {mustBeTextScalar} = 'wqficos';
        ficosOpt.iniFileName {mustBeTextScalar} = 'wqficos.ini';
        % file used as template for assigning parameters
        ficosOpt.settingsSourceName {mustBeTextScalar} = 'ficossettings.ini';
        % (optional) use this settings file directly 
        optFields.settingsFile {mustBeTextScalar} = '';
        % (optional) direct paths
        optFields.simFile {mustBeTextScalar} = '';
        optFields.iniFile {mustBeTextScalar} = '';
        optFields.settingsSource {mustBeTextScalar} = '';
        % check that files exist
        optFields.checkFiles {mustBeNumericOrLogical} = false;
        % include hdf5 options (paths to .hdf5 input files)
        optFields.hdf5Opt struct = struct;
    end

    if isempty(optFields.simFile)
        ficosOpt.simFile = fullfile(ficosOpt.simPath, ficosOpt.simFileName);
    else
        ficosOpt.simFile = optFields.simFile;
    end
    if isempty(optFields.iniFile)
        ficosOpt.iniFile = fullfile(ficosOpt.simPath, ficosOpt.iniFileName);
    else
        ficosOpt.iniFile = optFields.iniFile;
    end
    if isempty(optFields.settingsSource)
        ficosOpt.settingsSource = fullfile(ficosOpt.simPath,...
            ficosOpt.settingsSourceName);
    else
        ficosOpt.settingsSource = optFields.settingsSource;
    end
    if ~isempty(optFields.settingsFile)
        ficosOpt.settingsFile = optFields.settingsFile;
    end

    if ~isempty(optFields.hdf5Opt)
        % include also hdf5 options
        if isfield(optFields.hdf5Opt, 'hd_files')
            ficosOpt.hd_files = optFields.hdf5Opt.hd_files;
        end
        if isfield(optFields.hdf5Opt, 'loading')
            ficosOpt.loading = optFields.hdf5Opt.loading;
        end
    end

    if optFields.checkFiles
        % assert that files exist
        mustBeFile(ficosOpt.simFile);
        mustBeFile(ficosOpt.iniFile);
        if isfield(ficosOpt,'settingsFile')
            mustBeFile(ficosOpt.settingsFile);
        end
        mustBeFile(ficosOpt.settingsSource);
    end
end

