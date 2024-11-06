function [] = plotPostPred(f, t, wf_id, var, yPred, tbObs, options, figOpt)
%PLOTPOSTPRED Plot posterior predictive summary.
%   plotPostPred(f, t, wf_id, var) Plot posterior predictive summaries.
%
%   Required arguments:
%   f     : timetable with posterior predictions
%   t     : datetime column vector for the time period plotted
%   wf_id : codes for the water formations (e.g. intensive monitoring stations)
%           plotted
%   var   : variables plotted (e.g. 'DIN1', 'chla')
%
%   Optional arguments:
%   tbObs : (positional argument) time table of observations to be plotted
%           together with predictions. Default dataset used if argument is not
%           defined.
%   
%   Copyright (c) Karel Kaurila 2017 - 2024
%
arguments
    f (:,:) {mustBeNumeric}
    t (:,1) {mustBeA(t, 'datetime')}
    wf_id (:,1) {mustBeA(wf_id, 'categorical')}
    var (:,1) {mustBeA(var, 'categorical')}
    yPred (:,:) {mustBeNumeric} = NaN;
    tbObs (:,:) {mustBeA(tbObs, 'timetable')} = ...
        load_intens(true, 'old', 'all')
    options.Scale (1,:) {mustBeMember(options.Scale, {'linear', 'log'})} = 'linear';
    options.predObsAlpha (1,1) {mustBeNonnegative} = 0.1;
    % CI _edge_ color, separate from the face color inside the area
    options.intEdgeCol = 'blue';
    options.obsEdgeCol = 'cyan';
    % whether to plot median for posterior predictive observations 
    options.plotYPredMedian (1,1) {mustBeNumericOrLogical} = false;
    % for which variables to plot the CI for observations
    options.obsIntervalVars (1,:) {mustBeText} = {'A','FC','chla'};
    % Mark censored observations separately
    options.markCensored (1,1) {mustBeNumericOrLogical} = true;
    % Consider observations below simulator refugee values as censored
    options.ltRefAsCens (1,1) {mustBeNumericOrLogical} = true;
    % Order of plotted variables
    options.plotOrder (1,:) cell = ...
        {'DIN1','DIN2','DIP1','DIP2','chla','FC','A'};
    options.incTitle (1,1) {mustBeNumericOrLogical} = true;
    % median line width
    options.lwd (1,1) {mustBePositive} = 1.3;
    options.obsLwd (1,1) {mustBePositive} = 1.3;
    % observation marker sizes, non-censored and censored
    options.obsMrkSz (1,1) {mustBePositive} = 9;
    options.cenMrkSz (1,1) {mustBePositive} = 2;
    options.cenMrkColor char = 'black';
    % Export (saving) options
    options.exportFigs (1,1) {mustBeNumericOrLogical} = false;
    options.expOpt struct = struct('dir', 'figures/');
    % Axes options -------------
    options.minorTicks (1,1) {mustBeNumericOrLogical} = true;
    options.majorTicks ...
        {mustBeMember(options.majorTicks, {'yearly', 'quarterly'})} = 'yearly';
    options.yUnits (1,1) {mustBeNumericOrLogical} = false;
    options.yRot (1,1) {mustBeNumeric} = 0;
    options.compactYlab = false;
    options.axLwd = 0.5;
    options.majTkWidth (1,1) {mustBeNumeric} = 0;
    options.minTkWidth (1,1) {mustBeNumeric} = 0;
    options.gridLwd (1,1) {mustBeNumeric} = 0;
    options.minGridLwd (1,1) {mustBeNumeric} = 0;
    options.labFntSz = 10;
    options.tckFntSz = 8;
    options.alignYLabs (1,1) {mustBeNumericOrLogical} = false;
    options.manualYlabPos {mustBeNumeric} = [];
    % y tick label alignment
    options.YlabHorzAlign = 'center';
    % y tick label rotation
    options.YTkLabRot (1,1) {mustBeNumeric} = 0;
    % expand ylimits so until this much empty space at the top end
    % if 0, keep tickaligned yaxis
    options.YLimDelta (1,1) {mustBeNonnegative} = 0;

    % separate tile for ylab 
    options.ylabTile (1,1) {mustBeNumericOrLogical} = true;
    options.layoutOpt struct ...
        = struct('TileSpacing','tight','Padding','tight');
    options.innerLayoutOpt struct ...
        = struct();
    % when using separate tile for ylab, optionally also separate font size for
    % ylabel
    options.ylabFntSz (1,1) {mustBeNonnegative} = 0;
    % proportion of plot width to ylab width, in tiles 
    options.ylabTileProp (1,1) {mustBeInteger, mustBeNonnegative} = 14;
    options.useInnerTiles (1,1) {mustBeNumericOrLogical} = false;
    options.yTkLabLen (1,1) {mustBeNumeric} = 3;
    % put the major xtick label (year) between the ticks
    options.centerXtkLab (1,1) {mustBeNumericOrLogical} = false;
    % similarly, align y tick labels closer to the middle of axes
    options.hiddenYax (1,1) {mustBeNumericOrLogical} = false;

    % Figure options ---------------------
    figOpt.Units = 'pixels';
    figOpt.width = 2300;
    figOpt.height = 1300;

end


if(options.exportFigs)
    % Set exporting options
    expOpt = struct2opt(options.expOpt);
    exportFcn = @(fh, fname) exportFig(fh, 'fname', fname, expOpt{:});
end

wf_ids = unique(wf_id);
n_wf_ids = length(wf_ids);
plotVars = unique(var);

% Reordering plotVars:
typicalVars = options.plotOrder;
matchedTypVars = intersect(typicalVars,plotVars, 'stable');
otherVars = setdiff(plotVars, typicalVars);
orderedVars = [matchedTypVars(:); otherVars(:)];
plotVars = orderedVars;

nVars = length(plotVars);

% Include observed points if observation table given as an argument
includeObs = ~isempty(tbObs);

if includeObs
    % Only include observations within prediction period
    tbObs = cleanSimOutput(tbObs);
    tbObs = tbObs(tbObs.Time >= min(t) & tbObs.Time <= max(t), :);
    if ~ismember('variable', tbObs.Properties.VariableNames)
        tbObs = stackVars(tbObs, 'colName', 'observation', ...
            'idxName', 'variable');
    end

    if options.markCensored
        tbBounds = getCensureBounds('ltRefAsCens',options.ltRefAsCens);
        tbObs = innerjoin(tbObs, tbBounds, 'Keys', {'variable'}, ...
            'RightVariables', {'censorLB'});
        tbObs.censored = tbObs.observation < tbObs.censorLB;
    end

end

%% Quantiles

[f, fq, includePredQuantiles] = predQuantiles(f); 
% similarly for posterior predictive observations
[yPred, yqPred, includeYPredQuantiles] = predQuantiles(yPred);

%% Xticks (time)

% time vector used in plotting credible intervals
xt = unique(t, 'sorted');
xt_patch = [xt; flip(xt)];

yearsToPlot = unique(year(t), 'sorted');
if options.minorTicks
    % Monthly minor ticks and quarterly major ticks
    XTminor = arrayfun(@(x) datetime(x, 1:12,1), yearsToPlot, ...
        'UniformOutput',false);
    XTminor = [XTminor{:}];
end

switch options.majorTicks
    case 'quarterly'
        % Major ticks on 1. Jan, Apr, Jul and Oct
        XTmajor = arrayfun(@(x) datetime(x, [1, 4, 7, 10], 1), yearsToPlot, ...
            'UniformOutput', false);
        XTmajor = [XTmajor{:}];
    case 'yearly'
        % Major ticks on 1. Jan divide the years
        XTmajor = arrayfun( @(x) datetime(x, 1, 1), yearsToPlot);
end



if options.centerXtkLab
    % Year label in the middle of the year (1. July)
    % used in a workaround with separate axes:
    % https://se.mathworks.com/matlabcentral/answers/593956-how-to-write-in-the-middle-of-ticks
    XTmid = arrayfun( @(x) datetime(x, 7, 1), yearsToPlot);
end

%% Plotting options
layoutOpts = struct2opt(options.layoutOpt);
if isempty(options.innerLayoutOpt)
    innerLayoutOpt = layoutOpts;
else
    innerLayoutOpt = struct2opt(options.layoutOpt);
end

optMedian = {'Color', 'black', 'LineWidth', options.lwd};
optObsMedian = {'Color','magenta', 'LineStyle','--',...
    'LineWidth', options.obsLwd};

optInterval = {'blue', 'EdgeColor',options.intEdgeCol, ...
    'FaceColor','blue'};

optObsInterval = {'black', 'EdgeColor','none', ...
            'FaceColor', 'cyan', 'FaceAlpha', options.predObsAlpha};

optObsPoint = {'LineStyle','none','Marker','.','Color','red',...
    'MarkerSize', options.obsMrkSz};

% Censored observations not marked right now
optObsCensored = {'LineStyle','none',...
    'Marker','diamond',...
    'MarkerFaceColor',options.cenMrkColor,...
    'MarkerEdgeColor','none', 'MarkerSize', options.cenMrkSz};

legOpts = {'Location', 'northwest','Box','on'...
        'Orientation','horizontal'};

defFigPos = get(0, 'defaultFigurePosition');

figOpt2 = struct;
figOpt2.Units = figOpt.Units;
figOpt2.Position = [defFigPos(1:2), figOpt.width, figOpt.height];

figOptCell = struct2opt(figOpt2);


% Set tile structure
if options.ylabTile
    % Align Y axis labels by having them in a separate tile
    plotTileSpan = [1, options.ylabTileProp];
    ylabTileSpan = [1, 1];
else
    plotTileSpan = [1, 1];
    ylabTileSpan = [0, 0];
end

rowTileSpan = [1, plotTileSpan(2)+ylabTileSpan(2)];

%% Plotting

for i = 1:n_wf_ids
    figH = figure(figOptCell{:});
    tl = tiledlayout(nVars, rowTileSpan(2), layoutOpts{:});
    
    if options.useInnerTiles
        ylabLayout = tiledlayout(tl, nVars, ylabTileSpan(2), ... 
            innerLayoutOpt{:});
        ylabLayout.Layout.Tile = tilenum(tl, 1, 1);
        ylabLayout.Layout.TileSpan = [nVars, ylabTileSpan(2)];
        ylabParent = ylabLayout;
        
        innerLayout = tiledlayout(tl, nVars, rowTileSpan(2),...
            innerLayoutOpt{:});
        innerLayout.Layout.Tile = tilenum(tl, 1, 2);
        innerLayout.Layout.TileSpan = [nVars, rowTileSpan(2)];
        tileParent = innerLayout;
    else
        ylabParent = tl;
        tileParent = tl;
    end
    
    for j = 1:nVars
        
        

        if options.ylabTile
            ylabTile = nexttile(ylabParent, ylabTileSpan);
            axis(ylabTile, 'off');
        end
        
        ax = nexttile(tileParent, plotTileSpan);

        hiddenAxes = false;
        if options.hiddenYax || ...
                (options.centerXtkLab && j==nVars)
            % workaround for centered year labels
            % create an underlying copy axes with different xtick properties
            axUnder = cla(ax);
            if options.useInnerTiles
                ax = copyobj(axUnder, innerLayout);
            else
                ax = copyobj(axUnder, tl);
            end
            hiddenAxes = true;
        end
        indF = wf_id==wf_ids(i) & var==plotVars(j);

        plotPredictions(ax, f(indF), fq(indF,:), t(indF),...
            optInterval, optMedian);

        if(includeObs)
            indObs = tbObs.wf_id==wf_ids(i) & tbObs.variable==plotVars(j);
            if(ismember('censored', tbObs.Properties.VariableNames))
                indCen = indObs & tbObs.censored;
                plot(ax, tbObs(indCen,:),'observation', optObsCensored{:});
                
                indObs = indObs & ~tbObs.censored;
            end
            plot(ax, tbObs(indObs,:),'observation', optObsPoint{:});
        end
        
        YL = get(ax,'YLim');
        if ismember(plotVars(j), options.obsIntervalVars) || ...
                ismember('all', options.obsIntervalVars)
            plotPredictions(ax, yPred(indF), yqPred(indF,:), t(indF), ...
                optObsInterval, optObsMedian, ...
                plotMedian=options.plotYPredMedian);
        end


        if strcmp(options.Scale, 'linear')
            % Keep Ylimit as it was before posterior predictive for
            % observations
            ylim(ax, YL);
        end
        set(ax, 'YScale',options.Scale);
        
        % Setting up axes ticks  ----------
         
        % Major ticks
        set(ax, 'XTick', XTmajor);

        % Minor ticks
        if(options.minorTicks)
            ax.XAxis.MinorTickValues = XTminor;
            ax.XMinorGrid = 'on';
        end

        grid(ax, 'on');
        

        % Tick labels
        if(j < nVars)
            % Hide x-axis from all but the bottom plot
            ax.XAxis.Visible = 'off';
            set(ax, 'XtickLabel', []);
        else 
            if (options.minorTicks)
                ax.XAxis.MinorTick = 'on';
            end
        end

        if hiddenAxes && options.centerXtkLab
            % set tick labels on the secondary axes instead

            % invisible data point for datetimeruler to appear
            plot(axUnder, XTmajor(1), nan, 'x'); 
            set(ax, 'XTickLabel', []);
            if j==nVars
                set(axUnder,'XTick', XTmid, 'XTickLabel', year(XTmid), ...
                    'XLim', ax.XLim, 'YLim', ax.YLim);
            else
                set(axUnder, 'XtickLabel', []);
                axUnder.XAxis.Visible = 'off';
            end
            axUnder.YTick = [];
            axUnder.YLabel = [];
        end
    
        if hiddenAxes && options.hiddenYax
            % set ytick labels on the hidden y axis instead
            % to better control the label positioning
            ax.YAxis.Visible = 'off';
            % store original yticks and labels
            YTlab = ax.YTickLabel;
            YT = ax.YTick;
            set(ax, 'YTickLabel', []);
            charHeight = charHeightInDataUnits(ax, figH, options.tckFntSz);
            % move first tick label up by half a character
            YT(1) = YT(1)+0.5*charHeight;
            YLtmp = ax.YLim;
            if YT(end)+0.5*charHeight > YL(2) 
                % move last tick label down similarly
                YT(end) = YL(2) -0.5*charHeight;
            end
            set(axUnder, 'YTick', YT, 'YTickLabel', YTlab, ...
                'YLim', ax.YLim);


        end

        % axes linewidth - used for x(y)tick, edge line, box
        set(ax, 'LineWidth', options.axLwd);
        % Xtick width
        % Uses undocumented property:
        % https://se.mathworks.com/matlabcentral/answers/424549-how-could-i-change-the-linewidth-of-xtick-or-ytick-without-changing-the-box-or-edge-linewidth
        if options.majTkWidth > 0
            % if 0, use axLwd
            ax.XAxis.MajorTickChild.LineWidth = options.majTkWidth;
            % similarly for Ytick
            ax.YAxis.MajorTickChild.LineWidth = options.majTkWidth;
        end
        % set grid line widths
        if options.gridLwd > 0
            % if 0 use same as axes lwd
            ax.GridLineWidth = options.gridLwd;
        end
        if options.minGridLwd > 0
           ax.MinorGridLineWidth = options.minGridLwd;
        end

        
        xlabel(ax, '')
        
        % Y tick label rotation
        if options.YTkLabRot
            ax.YAxis.TickLabelRotationMode = 'manual';
            ax.YAxis.TickLabelRotation = options.YTkLabRot;
        end
        
        if options.YLimDelta
            expandYLim(ax, options.YLimDelta);
        end

        % Align ylabels
        ylab = varNameToYlab(string(plotVars(j)), options.yUnits, ...
                    'compact', options.compactYlab);
        ylabFcn = @(x) ylabel(x, ylab, 'Rotation', options.yRot);
        
        if options.ylabTile
            ylabTxt = text(ylabTile, 0, 0, ylab, 'Rotation', options.yRot);
             
            if all(options.ylabFntSz==0)
                ylabFntSz = options.labFntSz;
            else
                ylabFntSz = options.ylabFntSz;
            end
            ylabPos = [0.5, 0.5];
            if ~isempty(options.manualYlabPos) && ...
                    size(options.manualYlabPos,1) >= nVars && ...
                    size(options.manualYlabPos,2) == 2
                % set each ylabel position manually
                ylabPos = options.manualYlabPos(j, :);
            end
            ylabTxt.FontSize = ylabFntSz;
            set(ylabTxt, ...
                'Units', 'normalized', 'Position', ylabPos, ...
                'HorizontalAlignment', 'center');
            ylabel(ax, '');
        else
            if options.useInnerTiles
                ylabel(ax, '');
                ylabFcn(innerLayout);
            else
                ylabFcn(ax);
            end
    
            if options.alignYLabs && ~options.useInnerTiles
                % use same space for yticklabels so that ylabels align
                % e.g. %3g for 3 digits
                ytLabelFormat = sprintf('%%-%dg', options.yTkLabLen);
                ax.YAxis.TickLabelFormat = ytLabelFormat;
            end
        end

        % font sizes
        fontMult = options.labFntSz/options.tckFntSz;
        fontSzFcn = @(x) set(x, 'FontSize', options.tckFntSz, ...
            'LabelFontSizeMultiplier', fontMult);
        fontSzFcn(ax);
        if options.useInnerTiles
            innerLayout.YLabel.FontSize = options.tckFntSz;
        end
        if hiddenAxes
            % Match properties for hidden axes so that resizing doesn't break
            % things
            ylabel(axUnder, '');
            xlabel(axUnder, '');
            fontSzFcn(axUnder);
            set(axUnder, 'box', 'off');
            linkprop([ax, axUnder], {'Position', 'InnerPosition', ...
                'XLim', 'YLim'});
            % hide handle to the secondary axes
            axUnder.HandleVisibility = 'off';
        end
    end
    if(options.incTitle)
        title(tl, wf_id_to_title(wf_ids(i)))
    end

    %% Exporting
    if(options.exportFigs)
        fname = ['postPred_', wf_id_to_station(wf_ids(i))];
        exportFcn(gcf, fname);
    end

end

end
%% ---------------------------------

function [f, fq, incQ] = predQuantiles(f)
    incQ = false;
    fq = NaN(size(f,1), 2);
    if size(f, 2) == 3
        % lower quantile, median, upper quantile
        fq(:,1) = f(:,1);
        fq(:,2) = f(:,3);
        f = f(:,2);
    elseif size(f,2) == 2
        % Only lower and upper quantiles
        fq(:,1) = f(:,1);
        fq(:,2) = f(:,2);
        f = NaN(size(f,1), 1);
    end
end

function [] = plotPredictions(ax, f, fq, t, optInterval, optMedian, ...
    options)
    arguments
        ax (1,1) {mustBeA(ax, 'matlab.graphics.axis.Axes')} 
        f (:,1) {mustBeNumeric}
        fq (:,2) {mustBeNumeric}
        t (:,1) {mustBeA(t, 'datetime')}
        optInterval (1,:) cell
        optMedian (1,:) cell
        options.plotMedian {mustBeNumericOrLogical} = true;
    end
        y_patch = [fq(:,1); flip(fq(:,2))];
        t_patch = [t; flip(t)];
        patch(ax, t_patch, y_patch, optInterval{:});
        
        hold(ax, 'on');
        if options.plotMedian
             plot(ax, t, f, optMedian{:});
        end
end

function [] = expandYLim(ax, YLdelta)
    % expand YLim in an axis, when last ytick is close to upper ylimit
    % this way the top ytick label won't be so close to the next plot
    % and more space will be used for plots
    arguments
        ax {mustBeA(ax, 'matlab.graphics.axis.Axes')} = gca
        YLdelta (1,1) {mustBePositive} = 0.5;
    end
    %
    yTks = ax.YTick;
    YL = ax.YLim;
    if length(yTks) < 2
        return
    end
    tkGap = yTks(2)-yTks(1); 
    
    yTkRatio = (YL(2) - yTks(end))/tkGap;
    if yTkRatio < YLdelta
        YL(2) = yTks(end)+YLdelta*tkGap;
        ax.YLim = YL;
    end
end

function [charHeight] = charHeightInDataUnits(ax, fig, fontSz)
    % get the height of a character in the same units as the data in the current
    % axis
    arguments
        ax {mustBeA(ax, 'matlab.graphics.axis.Axes')} = gca;
        fig {mustBeA(fig, 'matlab.ui.Figure')} = gcf;
        fontSz {mustBePositive} = 7;
    end
    %
    figUnits = fig.Units;
    YL = ax.YLim;
    dataHeight = YL(2)-YL(1);
    if ~strcmp(ax.Units, 'normalized')
        error('Axes units must be normalied')
    end
    dataPerAxUnit = dataHeight/ax.Position(4);
    dataPerFigUnit = dataPerAxUnit/fig.Position(4);
    switch(figUnits)
        case 'centimeters'
            charHeightFigUnit = fontSz * 0.03528;
        case 'inches'
            charHeightFigUnit = fontSz/72;
        otherwise
            error('Only implented for cm or inches.')
    end
    charHeight = charHeightFigUnit * dataPerFigUnit;
end