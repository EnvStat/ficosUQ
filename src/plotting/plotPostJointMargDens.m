function [figH, exportPath] = plotPostJointMargDens(calRes, options, devOpt)
    %PLOTPOSTJOINTMARGDENS Plot parameter posterior (joint marginal) densities.
    %   plotPostJointMargDens(pathToCalibrationResults) Plots parameter
    %   posterior using the results in the given file.
    %
    %   [figH, exportPath] = plotPostJointMargDens(..., 'exportFig', true) 
    %   Saves the resulting figure to exportPath. figH is the figure handle to
    %   the plotted figure.
    %
    %   plotPostJointMargDens(..., 'passThrough', plotOpt) Set additional plotting
    %   options defined in structure plotOpt. See jointMargDens2D.m for what they 
    %   are.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        calRes {mustBeFile} = ...
            fullfile('results/calibration_example/calibration_results.mat');
        options.exportFig {mustBeNumericOrLogical} = false;
        % Additional plotting options passed through to the plotting functions
        options.passThrough struct = struct();
        % Developer options
        devOpt.trimTails (1,1) {mustBeNumericOrLogical} = false;
        devOpt.trimThresh (1,1) {mustBePositive} = log(2e3);
    end

    %% Load plotted samples from result file
    load(calRes, 'ws_rt', 'params');

    if devOpt.trimTails
        % Trim samples with low densities (high energy)
        load(calRes, 'energies');
        minEnergy = min(energies);
        iTrim = energies > minEnergy + devOpt.trimThresh;
        ws_rt = ws_rt(~iTrim, :);
    end

    plotArgs = namedargs2cell(options);
    [figH, exportPath] = jointMargDens2D(ws_rt, params, plotArgs{:});

    if ~isempty(exportPath)
        fprintf('Figure saved to %s\n', string(exportPath));
    end
end

