function [savePath] = getNewSaveFile(savepathOrig)
    %GETNEWSAVEFILE Generate a filepath that won't result overwriting existing
    %files.
    %   savePath = getNewSaveFile(savepathOrig) Checks if a file already exists
    %   at savepathOrig and if so returns a new path with timestamp attached.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        savepathOrig {mustBeTextScalar}
    end
    if isfile(savepathOrig)
            [parentDir, fName, fExt] = fileparts(savepathOrig);
            fName2 = sprintf('%s_%s%s', fName, getTimeStamp(), fExt);
            savePath = fullfile(parentDir, fName2);
    else
        savePath = savepathOrig;
    end
end

