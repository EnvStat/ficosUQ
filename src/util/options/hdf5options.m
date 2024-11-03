function [hdf5opt] = hdf5options(hdf5Dir, hdf5opt, options)
    %HDF5OPTIONS Paths to FICOS simulator's hdf5 input files.
    %   hdf5opt = hdf5options() Get default paths for the required hdf5 input
    %   files.
    %
    %   Named arguments should only be used by internal functions.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        % folder containing .hdf5 files
        hdf5Dir {mustBeTextScalar} = 'data/';
        % paths to .hdf5 hdf5 input files
        hdf5opt.loading {mustBeTextScalar} = ...
            fullfile(hdf5Dir, 'loading_Aurajoki.hdf5');
        hdf5opt.hd_files {mustBeTextScalar} = ...
            fullfile(hdf5Dir, 'hd_files.hdf5');
        % assert that files exist
        options.assertExist {mustBeNumericOrLogical} = false;
    end

    % ensure paths have the correct extension
    hdf5opt.loading = fileNameExt(hdf5opt.loading);
    hdf5opt.hd_files = fileNameExt(hdf5opt.hd_files);

    if options.assertExist
        % check that both files exist
        mustBeFile(hdf5opt.hd_files);
        mustBeFile(hdf5opt.loading);
    end
end

% --------
function [fpathOut] = fileNameExt(fpath, ext)
    % ensure given path has the specificed extension
    arguments
        fpath {mustBeText}
        ext {mustBeText} = '.hdf5';
    end
    [fpathRoot, fnameTmp] = fileparts(fpath);
    fpathOut = fullfile(fpathRoot, [fnameTmp, ext]);
end

