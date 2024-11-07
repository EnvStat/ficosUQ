function [pdfAx, cdfAx] = plotChlaDistSets(tbChlaSum, options)
    %PLOTCHLADISTSETS Plot density functions for mean seasonal chla.
    %   plotChlaDistSets(tbChlaSum) Plot pdf and cdf for seasonal chla using
    %   scenario summary table tbChlaSum.
    %
    %   See also SUMMARISECHLASCENARIOS RUNCATCHMENTSCENARIOS
    %
    %   Copyright (c) 2017-2024 Karel Kaurila
    %
    arguments
        tbChlaSum table
        % region for which summary is plotted
        options.region {mustBeTextScalar} = "Inner Archipelago";
        options.kernelWeights (1,1) {mustBeNumericOrLogical} = false;
        options.resetRng (1,1) {mustBeNumericOrLogical} = true;
        options.palette (:,3) {mustBeNumeric} = getDefaultPalette();
        % name of the scenario with no reductions
        options.BAUlegLabel {mustBeTextScalar} = 'Business As Usual';
        % use tiled layout
        options.useTiles {mustBeNumericOrLogical} = false;
        options.tlOpt struct = ...
            struct('TileSpacing', 'tight', 'Padding', 'none');
        options.figOpt struct = ...
            struct('Units', 'pixels', 'Position', [1920 700 650 283]);
        options.export (1,1) {mustBeNumericOrLogical} = true;
        options.expOpt struct = struct;
        options.legLoc {mustBeText} = {'northeast', 'southeast'};
        options.legFntSz (1,1) {mustBePositive} = 8;
        options.FntSz (1,1) {mustBePositive} = 10;
        options.TkLbFntSz (1,1) {mustBePositive} = 8;
        % Set axes labels again after plotting all individual lines
        % overrides labels set within plotScenarioChlaDist
        options.setAxesLabels (1,1) {mustBeNumericOrLogical} = false;
        options.includeUnits (1,1) {mustBeNumericOrLogical} = false;
    end

    if options.resetRng
        rng('default');
    end

    legLabels = [{options.BAUlegLabel}, ...
        compose('%d%% Load Reduction', 20:20:80)];

    weightOpt = {'UseWeightedKernel',options.kernelWeights};

    plotColors = options.palette;

    fontMult = options.FntSz/options.TkLbFntSz;
    axOpt = struct('FontSize', options.TkLbFntSz, ...
        'LabelFontSizeMultiplier', fontMult);

    axOpt = struct2opt(axOpt);
    legFontSize = options.legFntSz;
    xlimOpt = {};
    legOpt = {'Box', 'off', 'FontSize',legFontSize};

    legLocStr = string(options.legLoc);
    if isscalar(legLocStr)
        options.legLoc = [legLocStr; legLocStr];
    end

    if options.useTiles
        tlOpt = struct2opt(options.tlOpt);
        figOpt = struct2opt(options.figOpt);
        GESvalue = getGesValue(options.region);


        figPdf = figure(figOpt{:});
        tlPdf = tiledlayout(figPdf, 1, 1, tlOpt{:});
        pdfAx = nexttile(tlPdf);
        plotGESvalue(pdfAx, GESvalue);

        figCdf = figure(figOpt{:});
        tlCdf = tiledlayout(figCdf, 1, 1, tlOpt{:});
        cdfAx = nexttile(tlCdf);
        plotGESvalue(cdfAx, GESvalue);

        [pdfAx, cdfAx] = plotScenarioChlaDist(...
            tbChlaSum, options.region, ...
            'BAU', pdfAx, cdfAx,'LegendEntry', legLabels{1}, ...
            'Color', plotColors(1,:),...
            weightOpt{:});
    else
        [pdfAx, cdfAx] = plotScenarioChlaDist(tbChlaSum, options.region, ...
            'BAU', [], [],'LegendEntry', legLabels{1}, ...
            'Color', plotColors(1,:),...
            weightOpt{:});
    end


    catchLoads = 0.8:-0.2:0.2;
    for i = 2:length(legLabels)
        [pdfAx, cdfAx] = plotScenarioChlaDist(...
            tbChlaSum, options.region, ...
            'custom', pdfAx, cdfAx, 'CatchmentLoad', catchLoads(i-1),...
            'LegendEntry',legLabels{i}, 'Color', plotColors(i,:),...
            weightOpt{:});
    end

    if options.useTiles
        legend(pdfAx, 'location', options.legLoc{1}, legOpt{:});
        legend(cdfAx, 'location', options.legLoc{2}, legOpt{:});
    else
        legend(pdfAx, 'location', options.legLoc{1}, legOpt{:});
        legend(cdfAx,'location',options.legLoc{2}, legOpt{:});
        figPdf = pdfAx.Parent;
        figCdf = cdfAx.Parent;
    end
    xlim(pdfAx, xlimOpt{:});
    xlim(cdfAx, xlimOpt{:});
    set(pdfAx, axOpt{:});
    set(cdfAx, axOpt{:});

    if options.setAxesLabels
        setAxesLabels(pdfAx, 'pdf', 'enableUnits', options.includeUnits);
        setAxesLabels(cdfAx, 'cdf', 'enableUnits', options.includeUnits);
    end

    if(options.export)
        % export figures
        expOpt = options.expOpt;
        expOpt.fname = 'chlaPdf';
        expArgs = struct2opt(expOpt);
        exportFig(figPdf, expArgs{:});

        expOpt.fname = 'chlaCdf';
        expArgs = struct2opt(expOpt);
        exportFig(figCdf, expArgs{:});
    end
end
%% ---------------
function [] = setAxesLabels(ax, ylabType, options)
    % Set axes label strings
    arguments
        ax {mustBeA(ax, 'matlab.graphics.axis.Axes')} = gca;
        ylabType {mustBeTextScalar, ...
            mustBeMember(ylabType,{'pdf', 'cdf', ''})} = 'pdf';
        options.xlabStr {mustBeText} = "Summer chl\it{a}\rm";
        % Units for x label
        options.enableUnits {mustBeNumericOrLogical} = true;
        options.unitStr {mustBeTextScalar} = "{\mu}g chl\it{a}\rm/L";
        % ylab
    end
    % x label
    xlabStr = string(options.xlabStr);
    xlabStr = xlabStr(:);
    fullXlabStr = xlabStr;
    if options.enableUnits
        fullXlabStr = [xlabStr; string(options.unitStr)];
    end
    xlabel(ax, fullXlabStr);

    % y label
    setYlab = true;
    switch(ylabType)
        case 'pdf'
            ylabStr = sprintf("p(%s)", xlabStr);
        case 'cdf'
            ylabStr = sprintf("P(%s)", xlabStr);
        otherwise
            % only set x label
            setYlab = false;
    end

    if setYlab
        ylabel(ax, ylabStr);
    end

end
%% ------
function [] = plotGESvalue(ax, GESvalue)
    % plot the threshold for Good Environmental Status on the given axes
    arguments
        ax {mustBeA(ax, 'matlab.graphics.axis.Axes')} = gca;
        GESvalue = 3;
    end
    %
    xline(ax, GESvalue, 'LineStyle', '--', 'Color','black', ...
        'LineWidth', 1.5, 'HandleVisibility','off');
end
%% -------
function [GESvalue] = getGesValue(region)
    % returns the GES value for a given region
    arguments
        region {mustBeTextScalar} = "Inner Archipelago";
    end
    tb = getGesTable;
    GESvalue = tb.('GES value')(tb.region == string(region));
end
%% --------
function [pal] = getDefaultPalette()
    % Matlab's default palette
    % >> plot(1:3)
    % >> ax = gca;
    % >> ax.ColorOrder
    pal = [     0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
end