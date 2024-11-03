function [scOpt] = scenarioOptions(options)
    %SCENARIOMETADATA Define scenario options.
    %   scOpt = scenarioOptions() Get the definition of the default simulator
    %   scenario.
    %
    %   scOpt = scenarioOptions(...) Get specific scenario definition by setting
    %   named arguments.
    %
    %   Named arguments are
    %   'sampleID'     : indentifier for the sample used
    %   'scenarioID'   : identifier for the scenario
    %   'catchmentLoad : proportion for catchment area loadings relative to default
    %                    scenario
    %   'pointLoad'    : proportion for point source loadings relative to default
    %                    scenario
    %   'internalLoad' : proportion for internal loadings relative to default
    %                    scenario
    %   'atmosphericLoad' : proportion for atmospheric loadings relative to default
    %                    scenario
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.sampleID (1,1) {mustBeNumeric} = 1;
        options.scenarioID (1,1) {mustBeNumeric} = 1;
        options.catchmentLoad (1,1) {mustBeNumeric} = 1;
        options.pointLoad (1,1) {mustBeNumeric} = 1;
        options.internalLoad (1,1) {mustBeNumeric} = 1;
        options.atmosphericLoad (1,1) {mustBeNumeric} = 1;
    end
    scOpt = options;
end

