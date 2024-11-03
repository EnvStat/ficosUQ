function [scMetadata] = readScenarioMetadata(resultFile,options)
    %READSCENARIOMETADATA Read embedded scenario metadata from results.
    %   scMetadata = readScenarioMetadata(resultFile) Reads embedded scenario
    %   metadata from the given simulator output file.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        resultFile {mustBeFile};
        options.metadataName {mustBeTextScalar} = 'scenarioMetadata';
    end
    fileInfo = h5info(resultFile);
    dtsetNames = {fileInfo.Datasets.Name};
    if ~ismember(options.metadataName, dtsetNames)
        error('Scenario metadata  not found in %s, dataset %s missing.',...
            options.resultFile, options.metadataName);
    end
    encMtdt = h5read(resultFile, ['/',options.metadataName]);

    % format scenario information as struct array
    scMetadata = scenarioOptions;
    fnames = fieldnames(scMetadata);
    encMtdt = num2cell(encMtdt);
    encMtdt = encMtdt(:);
    scMetadata = cell2struct(encMtdt, fnames);
end

