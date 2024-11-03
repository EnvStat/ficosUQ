function [] = plotPostPredTimeseries(postPredResult,options)
    %PLOTPOSTPREDTIMESERIES Plot posterior predictive time series.
    %   plotPostPredTimeseries(postPredResult) Plots timeseries from posterior
    %   predictive results. postPredResult can either be timetable containing 
    %   posterior predictive summary obtained as the output to
    %   'postPredTimeseries', or it can be a file containing this timetable
    %   under a variable named 'tbPostPredSum'.
    %
    %   plotPostPredTimeseries() When called without arguments, plots figures
    %   for example results.
    %
    %   Additional (name, value) arguments:
    %   type - either 'figSet1' (default) or 'figSet2'. This determines the 
    %   the style of figures plotted.
    %
    %   save - true (default) / false, whether to save figures under 
    %   'figures/'.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        % which broad type of plot (sets multiple options)
        postPredResult  {mustBeA(postPredResult, {'timetable', 'char'})} = ...
            'results/calibration_example/postPredSummary.mat';
        options.type {mustBeMember(options.type, {'figSet1', 'figSet2'})} = ...
            'figSet1';
        options.save {mustBeNumericOrLogical} = true;
        options.showWarnings {mustBeNumericOrLogical} = false;
    end

    if ~isa(postPredResult, 'timetable')
        mustBeFile(postPredResult);
        load(postPredResult, 'tbPostPredSum');
        tbPPint = tbPostPredSum;
        clear('tbPostPredSum');
    else
        tbPPint = postPredResult;
        clear('postPredResult');
    end

    if ~options.showWarnings
        warning('off', 'MATLAB:table:ModifiedVarnamesUnstack');
    end

    if strcmp(options.type, 'figSet1')
        % Plot chla, A, FC, DIN1, DIP1
        plotOrder = {'chla', 'FC', 'A', 'DIN1', 'DIP1'};
        tbPPint(~ismember(tbPPint.variable, plotOrder),:) = [];

        % Plot years 2006 to 2009
        yearRange = [2006, 2009];
        tbPPint( year(tbPPint.Time) < yearRange(1), :) = [];
        tbPPint( year(tbPPint.Time) > yearRange(2), :) = [];
    end


    % Plotted observations
    tbObs = load_intens(true, 'old', 'summer');

    simCI = tbPPint{:,{'predLow','predMed','predHigh'}};
    timeVec = tbPPint.Time;
    wf_id_Vec = tbPPint.wf_id;
    varVec = tbPPint.variable;
    switch options.type
        case 'figSet1'
            % include post pred intervals for observations
            obsCI = tbPPint{:, {'obsLow', 'obsMed', 'obsHigh'}};
        case 'figSet2'
            % Do not include for this manuscript
            obsCI = tbPPint{:, {'obsMed'}};
    end

    % Shared plot options
    plotOpt = struct;
    plotOpt.Scale = 'linear';
    plotOpt.markCensored = true;
    plotOpt.ltRefAsCens = true;
    plotOpt.plotYPredMedian = false;
    plotOpt.predObsAlpha = 0.1;

    figOptS = struct();
    figOptS.Units = 'pixels';
    figOptS.width = 2300;
    figOptS.height = 1300;
    switch options.type
        case 'figSet1'
            plotOpt.obsIntervalVars = {'all'};
            plotOpt.plotOrder = plotOrder;
            plotOpt.yUnits = true;
            plotOpt.obsMrkSz = 9;
            plotOpt.cenMrkSz = 4;
            plotOpt.intEdgeCol = 'none';
            plotOpt.obsEdgeCol = 'none';
        case 'figSet2'
            plotOpt.obsIntervalVars = {};
            plotOpt.minorTicks = true;
            plotOpt.majorTicks = 'yearly';
            plotOpt.centerXTkLab = true;
            plotOpt.yUnits = true;
            plotOpt.yRot = 90;
            plotOpt.compactYlab = true;
            plotOpt.incTitle = false;
            % CI _edge_ color, separate from the face color inside the area
            plotOpt.intEdgeCol = 'blue';
            plotOpt.obsEdgeCol = 'none';
            % linewidths and marker sizes
            plotOpt.lwd = 0.6;
            plotOpt.obsLwd = 0.2;
            plotOpt.axLwd = 0.7;
            plotOpt.gridLwd = 1;
            plotOpt.minGridLwd = 0.5;
            plotOpt.obsMrkSz = 6;
            plotOpt.cenMrkSz = 1.5;
            plotOpt.labFntSz = 9;
            plotOpt.tckFntSz = 7;
            plotOpt.alignYLabs = true;
            plotOpt.ylabFntSz = 7;
            plotOpt.yTkLabLen = 3;
            % Adjust ylabel positions
            manYlPos = repmat([0.60, 0.40], 7, 1);
            manYlPos([1,2],2) = [0.2, 0.25]';
            manYlPos(5:7, 2) = [0.6, 0.6, 0.7]';
            plotOpt.manualYlabPos = manYlPos;

            % use hidden y axes to better position tick labels
            plotOpt.hiddenYax = true;
            % ytick label alignment

            plotOpt.useInnerTiles = true;
            plotOpt.layoutOpt = ...
                struct('TileSpacing','compact','Padding','tight');
            plotOpt.innerLayoutOpt = ...
                struct('TileSpacing','tight','Padding','tight');
            figOptS.Units = 'centimeters';
            figOptS.width = 17;
            figOptS.height = 9.35;
    end

    % export options
    if options.save
        plotOpt.exportFigs = true;
        exportDir = fullfile('figures', 'timeseries' ,getTimeStamp(2));
        mkdir(exportDir);
        plotOpt.expOpt = struct('dir', exportDir);
        plotOpt.expOpt.ext = {'png'};
    else
        plotOpt.exportFigs = false;
    end

    plotOpt = struct2opt(plotOpt);
    figOpt = struct2opt(figOptS);
    plotPostPred(simCI, timeVec, wf_id_Vec, varVec, obsCI, tbObs, ...
        plotOpt{:}, figOpt{:});
end