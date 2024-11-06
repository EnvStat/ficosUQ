function [figH, exportPath] = jointMargDens2D(xSamp, params, ...
        options, optExport, optExtra)
%JOINTMARGDENS2D Plot pairwise posterior joint marginal densities in a DxD 
% matrix off-diagonal, with (1D) marginals on the diagonal
%
% jointMargDens2D(samples, params) Plots posterior joint marginal densities from
% parameter samples into figures/jointPosterior_<timestamp> in both .png and
% .gif formats.
% 
% required arguments:
% samples - (n x d) matrix of posterior samples for the simulator parameters
% params  - structure containing calibration settings (parameters), needed for
% setting the names of the parameters etc.
% 
% optional (name, value) arguments:
% 'dir'  - path to folder where the resulting figure is saved (default:
% 'figures/')
% 'fname' - base filename used for the resulting figure (default:
% jointPosterior)
% 
% Other optional arguments are for fine-tuning the figure.
%
% Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    xSamp {mustBeNumeric}
    params struct
    options.paramUnits (1,1) {mustBeNumericOrLogical} = false;
    % Figure dimensions
    options.targetWidth (1,1) {mustBePositive} = 6.75; % ~= 17.20 cm
    options.widthToHeight (1,1) {mustBePositive} = 0.75;
    options.figUnits {mustBeMember(options.figUnits, ...
        {'pixels', 'inches', 'centimeters', 'normalized', 'points', 'characters'})} = 'inches';
    % Layout options -------
    options.ylabTile (1,1) {mustBeNumericOrLogical} = true;
    % proportion of ylabel width to total width
    options.ylabTileProp (1,1) {mustBeInteger, mustBePositive} = 18;
    options.mainSpacing ...
        {mustBeMember(options.mainSpacing, {'loose', 'compact','tight'})} = ...
        'tight';
    % Padding around the complete figure
    options.mainPadding ...
        {mustBeMember(options.mainPadding, {'loose', 'compact','tight','none'})} = ...
        'none';
    options.plotSpacing ...
        {mustBeMember(options.plotSpacing, {'loose', 'compact','tight'})} = ...
        'tight';
    options.plotPadding ...
        {mustBeMember(options.plotPadding, {'loose', 'compact','tight','none'})} = ...
        'tight';
    options.ylabSpacing ...
        {mustBeMember(options.ylabSpacing, {'loose', 'compact','tight'})} = ...
        'tight';
    options.ylabPadding ...
        {mustBeMember(options.ylabPadding, {'loose', 'compact','tight','none'})} = ...
        'tight';
    % Change colorbar position at the end
    options.moveCB (1,1) {mustBeNumericOrLogical} = true;
    % Axes options (font sizes) ------
    options.fontSize (1,1) {mustBePositive} = 10;
    options.tickLabSize (1,1) {mustBePositive} = 8;
    % Line options
    options.linewidth (1,1) {mustBePositive} = 1.5;
    % Reference dpi when converting units from inches to pixels manually
    options.useRefDpi (1,1) {mustBeNumericOrLogical} = false;
    options.referenceDpi (1,1) {mustBePositive} = 150;  
    % Export options
    options.exportFig (1,1) {mustBeNumericOrLogical} = false;
    optExport.fname {mustBeTextScalar} = 'jointPosterior';
    optExport.addTimeStamp {mustBeNumericOrLogical} = true;
    optExport.ext {mustBeText} = ...
        {'png','gif'};
    optExport.dir {mustBeTextScalar} = 'figures/';
    optExport.Resolution (1,1) {mustBePositive} = 600;
    % Pass through options from earlier function calls
    optExtra.passThrough struct = struct();
end

if ~isempty(optExtra.passThrough)
    options = assignFields(options, optExtra.passThrough);
    optExport = assignFields(optExport, optExtra.passThrough);
end
d = size(xSamp,2);
targetWidth = options.targetWidth;
targetHeight = options.widthToHeight*targetWidth;
pxPerInch = options.referenceDpi;
priorPlotN = 1000;
nBins = 30;

jointEdgeOpt = {'LineWidth',0.1,... % width only used if visible
                'EdgeColor','none'}; % remove these two to show edges
jointNorm = {'normalization', 'probability'}; % Joint histogram normalization method
lineOpt = {'lineWidth', options.linewidth, 'lineStyle', '-'};
postColor = {'Color', 'blue'};
priorColor = {'Color', 'black'};
xlimOpt = {'tight'};
tickLabelFontSize = options.tickLabSize;
labelFontSize = options.fontSize;
labMult = labelFontSize/tickLabelFontSize;
ylabOpt = {'VerticalAlignment', 'bottom'};

cbTitle = 'Posterior probability';
linkClim = true;

axOpt = {'FontSize',tickLabelFontSize, ...
    'LabelFontSizeMultiplier', labMult};

if(options.useRefDpi)
    % Older method for ensuring readable font sizes
    % Also creates figure on second monitor
    figOpt = {'Units','pixels',...
        'Position', [1920+700, 700, pxPerInch*targetWidth, pxPerInch*targetHeight]};
else
    figOpt = {'Units', options.figUnits, ...
        'Position', [0, 0, targetWidth, targetHeight]};
end

figH = figure(figOpt{:});

%% Set layout structure
if options.ylabTile
    % Additional tile for ylabel
    mainPlotTileSpan = [d, options.ylabTileProp-1];

    ylabTileSpan = [1, 1];
else
    ylabTileSpan = [0, 0];
    mainPlotTileSpan = [d, d];
end
totalTileSpan = [d, ylabTileSpan(2)+mainPlotTileSpan(2)];
% with nested layout, main plot can use the original indexing with
% [1,1] tile spans
plotTileSpan = [1, 1];

% Initialize axis limits as bounding box for samples
[lb, ub] = bounds(xSamp);
axisLims =  [lb', ub'];

paramLabels = cellfun(@(s)varNameToYlab(s, options.paramUnits), params.thetaNames,...
    'UniformOutput',false);

mainLayoutOpt = {'TileSpacing', options.mainSpacing, ...
    'Padding', options.mainPadding};

plotLayoutOpt = {'TileSpacing',options.plotSpacing,...
    'Padding',options.plotSpacing,...
    'TileIndexing','columnmajor'};
ylabLayoutOpt = {'TileSpacing', options.ylabSpacing,...
    'Padding', options.ylabPadding};


% main (outer) layout
tl = tiledlayout(figH, totalTileSpan(1), totalTileSpan(2), mainLayoutOpt{:});

if options.ylabTile
    % Additional column for y labels
    ylabCol = tiledlayout(tl, d-1, 1, ylabLayoutOpt{:});
    % Start from second row, since the top row only has a 1D marginal
    ylabCol.Layout.Tile = tilenum(tl, 2, 1);
    ylabCol.Layout.TileSpan = [d-1, 1];
    % Reserve space for ylabels first, 

    % inner layout for the rest of the plot
    axParent = tiledlayout(tl, d, d, ...
        mainLayoutOpt{:}, plotLayoutOpt{:});
    axParent.Layout.Tile = tilenum(tl, 1, 2);
    axParent.Layout.TileSpan = mainPlotTileSpan;
else
    set(tl, plotLayoutOpt{:});
    axParent = tl;
end


%% First plot the diagonal to update axis limits:
for i = 1:d
    tileInd = sub2ind([d d],i,i);
    
    ax = nexttile(axParent,tileInd, plotTileSpan);
    % Plot marginal density plot
    % Univariate kernel density estimation using 'kernel1', 
    % included in gpstuff
    [Pkernel, Xkernel] = kernel1(xSamp(:,i));
    axisLims(i,:) = [Xkernel(1), Xkernel(end)];
    xPrior = linspace(axisLims(i,1), axisLims(i,2), priorPlotN);
    pPrior = lognpdf(xPrior, params.mu_prior(i), params.sigma_prior(i));
    % Normalize pPrior, we are only interested in its shape
    pPrior = max(Pkernel)/max(pPrior)*pPrior;
    
    plot(xPrior, pPrior, priorColor{:}, lineOpt{:});
    hold on;
    plot(ax, Xkernel, Pkernel, postColor{:}, lineOpt{:});
    
    xlim(ax, xlimOpt{:});
    set(ax, axOpt{:});
    if(i==d)
        xlabel(ax, paramLabels{i});
    else
        xlabel(ax,'');
        xticklabels(ax, '');
    end
end
% Joint marginals
axJoint = zeros(nchoosek(d,2),1);
axInd = 0;
for i = 1:d
    for j = (i+1):d
        tileInd = sub2ind([d d],j,i);
        axInd = axInd+1;
        axJoint(axInd) = nexttile(axParent,tileInd, plotTileSpan);
        t1 = linspace(axisLims(i,1), axisLims(i,2), nBins);
        t2 = linspace(axisLims(j,1), axisLims(j,2), nBins);
        % Posterior joint marginal densities visualized through 
        % 2D histograms of the samples, normalized as probabilities
        histogram2(axJoint(axInd), xSamp(:,i), xSamp(:,j), ...
            t1, t2, ...
            'DisplayStyle', 'tile', ...
            jointNorm{:}, ...
            'FaceAlpha', 1,...
            'ShowEmptyBins','off', ...
            jointEdgeOpt{:});
        set(axJoint(axInd), 'XLim', axisLims(i,:), 'YLim', axisLims(j,:))
        colormap(axJoint(axInd), "hot");
        if(j==d)
            xlabel(axJoint(axInd), paramLabels{i});
        else
            xlabel(axJoint(axInd),'');
            xticklabels(axJoint(axInd), '');
        end
        if(i==1)
            if options.ylabTile
                % Hide original ylabel
                ylabel(axJoint(axInd), '');
                % Set ylabel on the ylabel tile as text element
                ylabTile = nexttile(ylabCol, tilenum(ylabCol, j-1, 1));
                axis(ylabTile, 'off');
                ylabText = text(ylabTile, 'String', paramLabels{j}, ...
                    'FontSize', labelFontSize, 'Rotation', 90, ...
                    'Units', 'normalized','Position', [0.5, 0.5], ...
                    'HorizontalAlignment', 'center');
            else
                ylabel(axJoint(axInd), paramLabels{j}, ylabOpt{:});
            end
        else
            ylabel(axJoint(axInd), '');
            yticklabels(axJoint(axInd),'');
        end
        set(axJoint(axInd), axOpt{:});
    end
end

if(linkClim)
    linkprop(axJoint,'CLim');
end



cb = colorbar(axJoint(end));

cbTileNum = tilenum(axParent, 1, d);
cbTileSpan = [d-1, 1];
% Create axes on the tile first, so that we can access it later without
% overwriting the colorbar
cbTile = nexttile(cbTileNum, cbTileSpan);
axis(cbTile, 'off');

cb.Layout.Tile = cbTileNum;
cb.Layout.TileSpan = cbTileSpan;
cb.Label.String = cbTitle;
cb.Label.FontSize = labelFontSize;

if options.moveCB
    % Move colorbar to from center of the column to the left
    cbPos = cb.Position;

    tilePos = cbTile.OuterPosition;
    cbPos(1) = tilePos(1)+tilePos(3)/4;
    cbPos(3) = tilePos(3)/5;

    cb.Location = 'manual';
    cb.Position = cbPos;
end

%% Set layout spacing options again
% otherwise right part of the figure might get clipped

if options.ylabTile
    set(tl, mainLayoutOpt{:});
    set(ylabCol, ylabLayoutOpt{:});
end
set(axParent, plotLayoutOpt{:});


%% Exporting
exportPath = '';
if(options.exportFig)
    if ~isfolder(optExport.dir)
        mkdir(optExport.dir)
    end
    exportArgs = struct2opt(optExport);
    exportPath = exportFig(figH, exportArgs{:});
end


end