function simOut = simOutToStruct(fpath,block, dtSets, sz)
    %SIMOUTTOSTRUCT Internal function for loading simulator outputs from hdf5
    %formatted results files.
    %
    % See 'readOutput'.
    %
    % Copyright (c) Karel Kaurila 2017 - 2024
    %
    readFcn = @(ds) h5read(fpath, sprintf('%s/%s', block, ds),...
        1, sz);
    simOut = struct();
    for i = 1:size(dtSets,1)
        simOut.(dtSets{i}) = readFcn(dtSets{i});
    end
    blkCat = categorical(extractAfter(string(block), '/Blocks/'));
    simOut.block = repmat(blkCat, sz, 1);

end