function [status] = change_loadings(varargin)
%CHANGE_LOADINGS Change nutrient loadings. Matlab wrapper for
%changeloadsnutrients.py.
% 
% change_loadings(source_file='data/loading.hdf5', ...
% int=1, catchment=1, point=1, atm=1) Creates a copy of the given nutrient 
% loadings file with loadings multiplied by given non-negative
% factors (loadings for which the factor equals 1 will not be changed).
%
% The arguments are:
%   'int'       : multiplies internal loadings by given factor
%   'catchment' : multiplies catchment area loadings by given factor
%   'point'     : multiplies point loadings by given factor
%   'atm',      : multiplies atmospheric loadings by given factor
%   'copy'      : either a path where the file with changed loadings will be,
%                 or logical [true] to use the default name: 
%                 ('<source_name>_copy.hdf5').
%
% For example, change_loadings(source_file='data/loading.hdf5',catchment=0.5)
% Creates a copy of 'data/loading.hdf5' to 'data/loading_copy.hdf5', in which 
% the catchment area nutrient loadings are halved.
%
%  Copyright (c) 2017 - 2024 Karel Kaurila
%

%% Initializing python arguments
load_types = {'intload', ... % Internal loading
              'vemala', ... % Catchment area loadings
              'atmdep', ... % Atmospheric loading
              'ah', ... % Point loadings from Ahvenanmaa
              'sm'}; % Point loadings from Archipelago Sea (Saaristomeri)
loads = ones(length(load_types),1);

python_stem = getPythonInterpeter();


%% Parsing inputs
ip = inputParser;
ip.CaseSensitive = false;

addOptional(ip, 'source_file', 'data/loading.hdf5', @isfile);

validFactors = @(x) (isscalar(x) && islogical(x) && all(~x) ) || ...
(~isscalar(x) && isvector(x) && length(x)==length(load_types) && all(x>=0));
addOptional(ip, 'factors', false, validFactors)

validCopy = @(x) isscalar(x) && ...
    (islogical(x) && all(x) ) ||  isstring(x) || ischar(x);

addParameter(ip, 'copy', true, validCopy);

validLoad = @(x) (isscalar(x)) && (isnumeric(x)) && (x>=0);
addParameter(ip, 'int', 1, validLoad);
addParameter(ip, 'catchment', 1, validLoad);
addParameter(ip, 'point', 1, validLoad);
addParameter(ip, 'atm', 1, validLoad);

% Which nutrient to change, or both
% 
addParameter(ip, 'nutrient', 'both', @(x)ismember(x, {'both','N','P'}));


% Path to python source folder
addParameter(ip, 'pythonDir', fullfile('src','python'), @isfolder);


parse(ip, varargin{:});

% path to python source
pythonDir = ip.Results.pythonDir;
pythonSrc = fullfile(pythonDir, 'changeloadsnutrients.py');


source_file = ip.Results.source_file;
copy_file = ip.Results.copy;
if(islogical(copy_file))
    if(copy_file)
        % Copy to <source_file>_copy.hdf5
        [srcDir, srcName] = fileparts(source_file);
        copyName = sprintf('%s_copy.hdf5', srcName);
        copy_file = fullfile(srcDir, copyName);
    else
        % Do not copy, change source file
        % Not allowed as input any more to prevent overwriting template files.
        copy_file = source_file;
    end
end

arg_factors = ip.Results.factors;
if(~isscalar(arg_factors))
    loads = arg_factors;
else
    loads(1) = ip.Results.int;
    loads(2) = ip.Results.catchment;
    loads(3) = ip.Results.atm;
    loads(4) = ip.Results.point;
    loads(5) = loads(4);
end

loadsN = loads;
loadsP = loads;
if strcmp(ip.Results.nutrient, 'N')
    % Only change nitrogen loads
    loadsP = ones(size(loads));
elseif strcmp(ip.Results.nutrient, 'P')
    % Only change phosphorus loads
    loadsN = ones(size(loads));
end

%% Set target_file
% Creat a copy of the source file
if(ischar(copy_file) || isstring(copy_file))
    target_file = copy_file;
    [cpSucc, cpMsg] = copyfile(source_file, target_file);
    if ~cpSucc
        warning(cpMsg)
    end
    if ~isfile(target_file) || ~endsWith(target_file, '.hdf5')
        error('change_loadings: Destination file %s must be a .hdf5 file!', target_file);
    end
end

if ~exist('target_file', 'var')
    error('change_loadings: Copy file path not defined.')
end

if strcmp(target_file, source_file)
    error('change_loadings: Attempting to overwrite loading source file: %s!',...
        source_file);
end

%% Changing loadings

for i = 1:length(load_types)
    if(loads(i)~=1)
        % changeloadsnutrients.py usage:
        % python3 changeloads.py <loadfile> <loadname> <N-multiplier> <P-multiplier>
        python_cmd = sprintf('%s %s %s %s %.9g %.9g', ...
            python_stem, pythonSrc, target_file, load_types{i}, loadsN(i), loadsP(i));
        [status, cmdout] = system(python_cmd);
        if(status)
            fprintf('change_loadings.m: Python command exited with status: %d\n',...
                status);
            error(cmdout);
            return
        end
    end
end


%% Everything ok
status = 0;
end

