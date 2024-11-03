function [] = writeScenarioMetadata(resultFile, scenarioDef, options)
    %WRITESCENARIOMETADATA Embed scenario metadata in a result.hdf5 file.
    %   writeScenarioMetadata(resultFile, scenarioDef) Embed metadata 
    %   describing scenario into the simulation result resultFile.
    %
    %   Use function scenarioOptions to construct the scenario definition 
    %   scenarioDef.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        resultFile {mustBeFile}
        scenarioDef struct = scenarioOptions();
        options.metadataName {mustBeTextScalar} = 'scenarioMetadata';
    end
    fileInfo = h5info(resultFile);
    dtsetNames = {fileInfo.Datasets.Name};
    % Scenario metadata encoded as numeric array, see scenarioOptions
    scMetadata = struct2array(scenarioDef);
    dtsetPath = ['/', options.metadataName];
    if ~ismember(options.metadataName, dtsetNames)
        % No dataset for metadata created yet, creating it first
        h5create(resultFile, dtsetPath, size(scMetadata));
    end
    h5write(resultFile, dtsetPath, scMetadata);
end

