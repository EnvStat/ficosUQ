function [status,md5sum] = setupLoadingFile(options)
%SETUPLOADINGFILE Initializes nutrient loadings file.
%   setupLoadingFile() Checks that the nutrient loadings file exists and
%   verifies its integrity by checking its md5sum.
%
%   Used internally when initializing FICOS calibration.
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    options.hdf5path {mustBeTextScalar} = 'data/';
    options.filename {mustBeTextScalar} = 'loading';
    options.fullpath {mustBeTextScalar} = '';
end

validNames = {'loading', ... % Default loading file
              'loading_Aurajoki'}; % Loadings specific to Aurajoki estuary

status = -1; % not yet done
fullpath = '';

%% Input parsing ---

% Using full path if given, otherwise constructing it from individual parts
if ~isempty(options.fullpath)
    fullpath = options.fullpath;
    [hdf5path, filename, fExt] = fileparts(fullpath);
    if isempty(fExt)
        fullpath = [fullpath, '.hdf5'];
    elseif ~strcmp(fExt, '.hdf5')
        error('%s is not an .hdf5 file.\n', fullpath);
    end
else
    hdf5path = options.hdf5path;
    [~, filename, fExt] = fileparts(options.filename);
    if ~isempty(fExt) && ~strcmp(fExt, '.hdf5')
        warning('File %s is not a .hdf5 file', filename);
    end
    fullpath = fullfile(hdf5path, ...
                        [filename '.hdf5']);
end

%% 

if(isfile(fullpath))
    fprintf('File %s.hdf5 already exists, checking md5sum.\n', filename)
    md5sum = getMd5sum(fullpath);
    validMd5sum = validateLoadingFile(filename, md5sum);
    status = ~validMd5sum;
    if(~validMd5sum)
        error('%s.hdf5 - Invalid md5sum: %s.', filename, char(md5sum));
    end
elseif strcmp(filename, 'loading')
    %status = 1;
    % Base loading.hdf5 not found - cannot proceed further
    error('Base loading file %s not found under %s.', filename, hdf5path);
else
    source_file = fullfile(hdf5path, 'loading.hdf5');
    fprintf('Base loading file: %s\n', source_file);
    if ismember(filename, validNames)
        fprintf('Creating predefined loading file %s from base loading file\n',...
            filename);
        if strcmp(filename, 'loading_Aurajoki')
            % internal nitrogen loads reduced to 0.3 of base scenario
            change_loadings(source_file, 'copy', fullpath, ...
                            'nutrient', 'N', 'int', 0.3);
        end
    else
        % Create a copy of base file with given name
        [copySuccess, copyMsg, copyMsgId] = copyfile(source_file, fullpath);
        if(~copySuccess)
            warning(copyMsg);
            status = copyMsgId;
        else
            status = 0;
        end
    end
end

% Everything ok
status = 0;

if(nargout>=2)
    % Return md5sum of the loading file
    md5sum = getMd5sum(fullpath);
end


end
%% -----------------------------------------------------
function [out] = validateLoadingFile(fname, md5sum)
% Check that given loading file is correct with md5sum

% md5sums for known files
switch fname
    case 'loading' % base loading file
        out = strcmp(md5sum, 'd749d559aa16619757255394a6abe2d1');
    case 'loading_Aurajoki'
        out = strcmp(md5sum, '2b27168bc1b9d3d3cfd929463f2b2f63');
    otherwise
        warning(['Loading file %s not specific to any scenario. ' ...
            'It has md5sum %s.'], ...
            fname, char(md5sum));
        out = 2;
end
if(out)
    fprintf('Md5sum ok.\n');
end

end

function [md5sum, status, out] = getMd5sum(path)
cmd = ['md5sum ' path];
[status, out] = system(cmd);
if(~status)
    tmp = strsplit(out, ' ');
    md5sum = tmp{1};
else
    error('Could not check md5sum for %s\n Status: %d\n Output: %s\n',...
        path, status, string(out));
end


end