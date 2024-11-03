function [simOutput] = run_sim(theta, thetaNames, ...
        ind, wid, options, outputOpt, overrideOpt)
% RUN_SIM Launch a single FICOS simulation.
%   [simOut] = run_sim() Runs FICOS simulator with default parametrisation and
%   returns the predictions in a timetable.
%
%   [simOut] = run_sim(theta, thetaNames) Runs FICOS with parametrisation theta,
%   with parameters named thetaNames.
%
%   Additional settings can be defined through named arguments. The preferred
%   way of setting these is to use the helper function 'ficosOptions' as
%   follows:
%   simOut = run_sim(theta, thetaNames, 'ficosOpt', ficosOptions(...)).
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    theta             (1,:) {mustBeNumeric} = [10 15 10 0.1353 0.1353];
    thetaNames        (1,:) cell = ...
        {'Klight', 'LightN2fix', 'LightThres', 'RAmax', 'RFCmax'};
    % arguments ind and wid only meant for internal use
    ind               (1,1) {mustBeNumeric} = 1;
    wid               (:,:) = [];
    options.hd_files  (1,:) {mustBeTextScalar} = 'data/hd_files.hdf5';
    options.loading   (1,:) {mustBeTextScalar} = 'data/loading.hdf5';
    options.local_dir (1,:) {mustBeTextScalar} = '.';
    options.simPath   (1,:) {mustBeTextScalar} = 'bin/ficos/';
    options.debugPrints (1,1) {mustBeNumericOrLogical} = false;
    options.settingsSource (1,:) {mustBeTextScalar} = 'bin/ficos/ficossettings.ini';
    % settings file with parameters already set, which
    % overrides arguments theta, thetaNames and settingsSource
    options.settingsFile {mustBeTextScalar} = '';
    options.iniFile (1,:) {mustBeTextScalar} = 'bin/ficos/wqficos.ini';
    % output options
    outputOpt.outputRoot (1,:) {mustBeTextScalar} = '';
    outputOpt.outputPath (1,:) {mustBeTextScalar} = '';
    outputOpt.rmOldOutput (1,1) {mustBeNumericOrLogical} = false;
    outputOpt.outputFcn {mustBeA(outputOpt.outputFcn, 'function_handle')} = ...
        @(x) readOutput(x);
    % alternatively define simulator paths through a struct array
    % (overrides paths set above)
    overrideOpt.ficosOpt struct = struct;
end

%% Input parsing

% override path settings from ficosOpt, if given
if ~isempty(overrideOpt.ficosOpt)
    options = assignFields(options, overrideOpt.ficosOpt);
end

% Handling default arguments for output:
if isempty(outputOpt.outputRoot)
    % Default outputRoot - main folder for outputs
    % contains subfolders for multiple simulations
    outputOpt.outputRoot = fullfile('output/');
    if ~isempty(wid)
        outputOpt.outputRoot = fullfile(outputOpt.outputRoot, ...
            num2str(wid));
    end
end

if isempty(outputOpt.outputPath)
    % Default outputPath - contains outputs for a single simulation
    outputOpt.outputPath = fullfile(outputOpt.outputRoot, ...
        num2str(ind));
end

if outputOpt.rmOldOutput
    rmdir(outputOpt.outputPath, 's');
    pause(0.1);
end

% Suppress warning
warning('off', 'MATLAB:MKDIR:DirectoryExists');

mkdirSuccess = mkdir(outputOpt.outputPath);
if ~mkdirSuccess
    error('Could not create folder %s\n', outputOpt.outputPath);
end

if(options.debugPrints)
    disp(hd_files)
    disp(loading)
end

if isempty(options.settingsFile)
    settingsPath = fullfile(outputOpt.outputPath, ...
        sprintf('ficossettings%d.ini',ind));
    writeParamsToSettings(settingsPath, theta, thetaNames, ...
        options.settingsSource);
else
    settingsPath = options.settingsFile;
end
logOutFile = fullfile(outputOpt.outputPath, 'simOut.log');
logErrFile = fullfile(outputOpt.outputPath, 'simErr.log');

%% Compose simulator command
simCmd = fullfile(options.simPath, 'wqficos');
if(isempty(wid))
    cmd = sprintf('%s %s %s %s %s %s test 0 0 1>%s 2>%s',...
    simCmd, settingsPath, outputOpt.outputPath, ...
    options.hd_files, options.loading, options.iniFile, ...
    logOutFile, logErrFile);
else
    cmd = sprintf('%s %s %s input/w%d/%s input/w%d/%s wqficos.ini test 0 0 1>logs/sim_out%d.txt 2>logs/sim_err%d.txt',...
    simCmd, settingsPath, outputOpt.outputPath, ind, options.hd_files, ...
    ind, options.loading, wid ,wid);
end
[status, ~] = system(cmd);
if(status~=0)
   fprintf('Error with cmd:\n %s\n', cmd)
   error('Simulation ended with status %s.',string(status))
end
pause(1);

if nargout>= 1
    simOutput = outputOpt.outputFcn(outputOpt.outputPath);
end

end
