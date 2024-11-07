function [pdfAx, cdfAx, figECDF] = plotScenarioChlaDist(varargin)
    %PLOTSCENARIOCHLADIST Plot chlorophyll-alpha distribution for a given
    %scenario and region / water formation



    %% Input parsing

    % Default values and sets of valid values
    defaultResultTablePath = ...
        'results/scenarios/loads_int_catch_point/16-Sep-01:43:55/summaries_20231005.mat';
    defaultRegion = 'Inner Archipelago';
    defaultScenario = 'custom';

    defaultInternalLoad = 1;
    defaultCatchmentLoad = 1;
    defaultPointLoad = 1;

    % Valid scenarios:
    % Scenario others than custom override specified loadings
    % custom                        : specify loadings, if custom selected with
    %                                 no loadings specificed, same as BAU
    % BAU (Business as usual)       : no reduction in loadings
    % BSAP(Baltic Sea Action Plan)  : TODO - find corresponding loadings
    validScenarios = {'BAU','custom', 'BSAP'};

    % Validator functions ------------------------

    validResultTable = @(x) (istable(x)) || (isfile(x));
    validRegion = @(s) (isstring(s) || ischar(s) || iscategorical(s));
    fValidScenario = @(s) (isstring(s) || ischar(s)) && ...
        ismember(s,validScenarios);
    validLoad = @(x) isscalar(x) && (  (iscategorical(x)) || ...
        (isnumeric(x) &&  x>= 0 && x <= 1));
    validLoadVec = @(x) isvector(x) && ...
        (size(x,1) == 3 || size(x,1) == 4) && ...
        iscategorical(x) || ( isnumeric(x) && all(x>=0) && all(x<=1));
    validAx = @(x) isempty(x) || isa(x,'matlab.graphics.axis.Axes');
    % Optional arguments -----------------------
    ip = inputParser;
    ip.addOptional('resultTable', defaultResultTablePath, validResultTable);
    ip.addOptional('region', defaultRegion, validRegion);
    ip.addOptional('scenario', 'custom', fValidScenario);
    ip.addOptional('pdfAx', [], validAx);
    ip.addOptional('cdfAx', [], validAx);

    % Parameters --------------------------
    % Individually, or with one vector - vector overrides others
    defaultInternalLoad = 1;
    defaultCatchmentLoad = 1;
    defaultPointLoad = 1;
    % Atmospheric loads included here already, but as of 2023 Oct 5, we have
    % fixed this at 1, i.e. no reduction - also not included in result table
    % yet
    includeAtmLoad = false;
    defaultAtmosphericLoad = 1;
    defaultLoadVec = [1 1 1 1]; % Same order as above

    ip.addParameter('InternalLoad', defaultInternalLoad, validLoad);
    ip.addParameter('CatchmentLoad', defaultCatchmentLoad, validLoad);
    ip.addParameter('PointLoad', defaultPointLoad, validLoad);
    ip.addParameter('AtmosphericLoad', defaultAtmosphericLoad, validLoad);
    ip.addParameter('LoadVector', defaultLoadVec, validLoadVec);
    ip.addParameter('LegendEntry', '', @(x) ischar(x) || istring(x) );
    % Older 'Color' argument
    ip.addParameter('gesColor', []);
    % New 'Color' argument used for current function line
    ip.addParameter('Color', []);

    ip.addParameter('UseWeightedKernel', false)

    % density estimation method, either 'kernel' for kernel density estimation
    % or 'icdf' for inverse cdf sampling
    ip.addParameter('densEstMethod', 'kernel');

    loadParams = {'InternalLoad', 'CatchmentLoad', 'PointLoad', ...
        'AtmosphericLoad', 'LoadVector'};

    % Parsing --------------------------
    parse(ip, varargin{:});

    useWeightedKernel = ip.Results.UseWeightedKernel;
    densEstMethod = ip.Results.densEstMethod;

    region = ip.Results.region;
    if(~iscategorical(region))
        if(ischar(region))
            region = cellstr(region);
        end
        region = categorical(region);
    end
    scenario = ip.Results.scenario;

    T_summaries = table;
    resultTable = ip.Results.resultTable;

    gesColor = ip.Results.gesColor;
    optColor = struct;
    if(~isempty(gesColor))
        optColor.GEScolor = gesColor;
    end

    if ~isempty(ip.Results.Color)
        optColor.plotColor = ip.Results.Color;
    elseif ~isempty(gesColor)
        optColor.plotColor = gesColor;
    end

    optColor = struct2opt(optColor);

    if(istable(resultTable))
        T_summaries = resultTable;
        clear('resultTable');
    elseif(isfile(resultTable))
        load(resultTable, 'T_summaries');
    end

    % Append plots to existing axes or create new ones


    [pdfAx,pdfPlotNew] = checkPlotAx(ip.Results.pdfAx);
    [cdfAx,cdfPlotNew] = checkPlotAx(ip.Results.cdfAx);

    legendEntry = ip.Results.LegendEntry;

    % Parsing load parameters
    if(ismember('Atmospheric load', T_summaries.Properties.VariableNames))
        includeAtmLoad = true;
    end
    loadColumns = {'Internal load','Catchment area load','Point load'};
    if(includeAtmLoad)
        loadColumns = [loadColumns(:)', {'Atmospheric load'}];
    end

    % Collect loadings that appear in data
    loadCategories = categories(T_summaries{:,loadColumns});

    loadVector = categorical(ones(size(loadColumns)));
    if(strcmp(scenario, 'custom') && ...
            any(~ismember(loadParams, ip.UsingDefaults)))
        if(~ismember('LoadVector',ip.UsingDefaults))
            loadVector = ip.Results.LoadVector;
            if(~iscategorical(loadVector))
                loadVector = categorical(loadVector);
            end
        else
            loadVector(1) = categorical(ip.Results.InternalLoad);
            loadVector(2) = categorical(ip.Results.CatchmentLoad);
            loadVector(3) = categorical(ip.Results.PointLoad);
            if(includeAtmLoad)
                loadVector(4) = categorical(ip.Results.AtmosphericLoad);
            end
        end
    end

    loadVector = setcats(loadVector, loadCategories);
    if(any(isundefined(loadVector)))
        error('Specified scenario loads do not appear in result table.');
    end

    %% Filter to specified scenario
    rf = rowfilter(T_summaries);
    rf = (rf.region == region);
    for li = 1:length(loadColumns)
        rf = rf & (rf.(loadColumns{li}) == loadVector(li));
    end

    chlaCDFavailable = true;
    if ~ismember('chla CDF', T_summaries.Properties.VariableNames)
        % no weight information, using uniform weights
        chlaCDFavailable = false;
        useWeightedKernel = false;
    end

    if(chlaCDFavailable)
        T_plot = T_summaries(rf, {'GES value', 'Pr(GES)','chla CDF'});
        T_CDF = T_plot.('chla CDF'){1};
        GES_value = T_plot.('GES value');

        chlaVec = T_CDF{:,1};
        if useWeightedKernel
            % Convert cdf back to weights
            weightVec = T_CDF{:,2};
            weightVec(2:end) = weightVec(2:end)-weightVec(1:(end-1));
        else
            weightVec = ones(size(T_CDF{:,2}));
        end
        T_plot = array2table([chlaVec, weightVec],'VariableNames',...
            {'chla', 'weight'});
    else
        T_plot = T_summaries(rf, {'chla','GES value'});
        T_plot.weight = ones(size(T_plot.chla));
        GES_value = T_plot{1,'GES value'};
    end



    if(cdfPlotNew)
        xline(cdfAx, GES_value, 'LineStyle', '--', ...
            'Color','black', 'LineWidth', 1.5, ...
            'HandleVisibility','off');
    end
    if(pdfPlotNew)
        xline(pdfAx, GES_value, 'LineStyle', '--', ...
            'Color','black', 'LineWidth', 1.5, ...
            'HandleVisibility','off');
    end
    %% Sample from CDF
    if strcmp(densEstMethod, 'icdf')
        cdfVec = T_CDF.CDF;
        chlaVec = T_CDF.chla;
        nCdfVec = length(cdfVec);
        N_PDFSAMP = 10*nCdfVec;

        % Inverse CDF sampling
        chlaSamp = invCDFsamp(chlaVec, cdfVec, N_PDFSAMP);
        [pKernel, xKernel] = kernel1(chlaSamp);
    else
        [pKernel, xKernel] = kernelDensFromWeights(T_plot.('chla'), T_plot.('weight'));
        pKernel = reshape(pKernel, length(pKernel),1);
        xKernel = reshape(xKernel, length(xKernel),1);
    end


    pdfAx = plotPDF(xKernel, pKernel, pdfAx, GES_value,...
        'LegendEntry', legendEntry, optColor{:});


    %% Plot CDF


    if(nargout>=3)
        figECDF = nFig();
        axECDF = axes(figECDF);
        axECDF = plotCDF(T_CDF.chla, T_CDF.CDF, axECDF, GES_value, ...
            'stairs', true, optColor{:});
    end

    % Smoothed ECDF using kernel densities
    cdfKernel = cumsum(pKernel)/sum(pKernel);

    cdfAx = plotCDF(xKernel, cdfKernel, cdfAx, GES_value, ...
        'stairs', false, 'LegendEntry', legendEntry,...
        optColor{:});

end
% --------------------------------------------------------------
function [xSamp] = invCDFsamp(xVec, cdfVec, nSamp)
    uCDF = rand(1,nSamp);
    [II, JJ] = find( uCDF <= cdfVec);
    [~, IJ] = unique(JJ, 'first');
    cdfInds = II(IJ);
    upperInds = min(cdfInds, size(cdfVec,1));
    lowerInds = max(cdfInds-1, 1);
    xUb = xVec(upperInds);
    xLb = xVec(lowerInds);
    xSamp = xLb+(xUb-xLb).*rand(nSamp,1);
end
function [ax] = plotCDF(xVec, cdfVec, varargin)
    ip = inputParser;
    validGES = @(x) isempty(x)|| (isscalar(x) && isnumeric(x) && x>=0);
    validAx = @(x) isempty(x) || isa(x,'matlab.graphics.axis.Axes');
    validColor = @(x) isempty(x) || ischar(x) || isstring(x) || ...
        (isnumeric(x) && isvector(x) && (length(x)==3 || length(x)==4));

    ip.addOptional('ax', gca, validAx);
    ip.addOptional('GES', [], validGES);


    ip.addParameter('stairs', false, @(x)all(x)||all(~x));
    ip.addParameter('GEScolor', [], validColor);
    ip.addParameter('plotColor', [], validColor);
    ip.addParameter('LegendEntry', '');

    parse(ip,varargin{:});

    GES = ip.Results.GES;
    GEScolor = ip.Results.GEScolor;
    plotColor = ip.Results.plotColor;
    legendEntry = ip.Results.LegendEntry;

    newPlot = false;
    if(isempty(ip.Results.ax))
        figH = nFig();
        ax = gca();
        newPlot = true;
        hold on;
    else
        ax = ip.Results.ax;
        axes(ax);
        hold on;
    end


    if(ax.ColorOrderIndex==1)
        ax.ColorOrder = turbo(20);
        newPlot = true;
    end
    if isempty(GEScolor)
        GEScolor = ax.ColorOrder(ax.ColorOrderIndex,:);
    end
    if isempty(plotColor)
        plotColor = ax.ColorOrder(ax.ColorOrderIndex,:);
    end

    smoothPlot = ~all(ip.Results.stairs);
    if(~smoothPlot)
        [yy, xx] = stairs(cdfVec, xVec);
    else
        xx = xVec;
        yy = cdfVec;
    end
    plot(xx, yy, 'lineStyle', '-', 'LineWidth', 1.5,...
        'DisplayName',legendEntry, 'Color', plotColor);

    XL = ax.XLim;
    YL = ax.YLim;

    XL(1) = min(XL(1), xVec(1));
    XL(2) = max(XL(2), xVec(end));
    xlim(ax,XL);

    xMargin = 0.6;
    yMargin = 0.1;


    if(~isempty(GES) && GES > xVec(1))
        iGES = find(xVec>=GES, 1, 'first');
        if(~isempty(iGES) && iGES>1 && smoothPlot)
            % Linear interpolation
            dCdf = cdfVec(iGES)-cdfVec(iGES-1);
            dX = xVec(iGES)-xVec(iGES-1);
            PrGES = cdfVec(iGES-1) + (GES-xVec(iGES-1))*dCdf/dX;
        elseif(isempty(iGES))
            iGES = length(xVec);
            PrGES = 1;
        else
            PrGES = cdfVec(iGES);
        end

        if(PrGES >= 0.5)
            yGES = PrGES-yMargin;
            xGES = GES-1*xMargin;
        else
            yGES = PrGES+yMargin;
            xGES = GES-1*xMargin;
        end
        plot([-5, GES], PrGES*ones(1,2), 'Color', GEScolor, ...
            'LineStyle','-.','HandleVisibility','off');
        PrStr = sprintf('Pr(GES) = %3.2g', PrGES);
        text(ax, xGES, yGES, PrStr, 'Color','black');
    end
    XL(1) = min(min(XL(1), GES-1*xMargin), xVec(1));
    XL(2) = max(max(XL(2), GES+1*xMargin), xVec(end));
    xlim(ax,XL);
    ylim(ax, [0, 1])
    xlabel('Chlorophyll \it{a}\rm')
    ylabel('P(Chlorophyll \it{a}\rm)')
end
function [ax] = plotPDF(xVec, pVec, varargin)
    ip = inputParser;
    validGES = @(x) isempty(x)|| (isscalar(x) && isnumeric(x) && x>=0);
    validAx = @(x) isempty(x) || isa(x,'matlab.graphics.axis.Axes');
    validColor = @(x) isempty(x) || ischar(x) || isstring(x) || ...
        (isnumeric(x) && isvector(x) && (length(x)==3 || length(x)==4));

    ip.addOptional('ax', gca, validAx);
    ip.addOptional('GES', [], validGES);

    ip.addParameter('GEScolor', [], validColor);
    ip.addParameter('plotColor', [], validColor);
    ip.addParameter('LegendEntry', '');
    parse(ip,varargin{:});

    GES = ip.Results.GES;
    legendEntry = ip.Results.LegendEntry;
    newPlot = false;

    if(isempty(ip.Results.ax))
        figH = nFig();
        ax = gca();
        newPlot = true;
    else
        ax = ip.Results.ax;
        axes(ax);
    end

    GEScolor = ip.Results.GEScolor;
    plotColor = ip.Results.plotColor;

    if(ax.ColorOrderIndex==1)
        if isMATLABReleaseOlderThan("R2023b")
            ax.ColorOrder = turbo(20);
        else
            colororder(ax, "gem12");
        end
        newPlot = true;
    end
    hold(ax, 'on');
    if isempty(GEScolor)
        GEScolor = ax.ColorOrder(ax.ColorOrderIndex,:);
    end
    if isempty(GEScolor)
        plotColor = ax.ColorOrder(ax.ColorOrderIndex,:);
    end
    plot(ax, xVec, pVec, 'LineStyle','-','LineWidth',1.5, ...
        'DisplayName',legendEntry, 'Color', plotColor);

    XL = ax.XLim;
    YL = ax.YLim;
    xMargin = range(XL)*0.05;
    if(~isempty(GES))
        if(GES > xVec(1))
            iGES = find(GES <= xVec,1,'first');
            if( isempty(iGES) || iGES==length(pVec))
                iGES = length(pVec);
                yGES = pVec(end);
            else
                dP = pVec(iGES)-pVec(iGES-1);
                dX = xVec(iGES)-xVec(iGES-1);
                yGES = pVec(iGES-1) + (GES-xVec(iGES-1))*dP/dX;
            end
            xpatch = [xVec(1:(iGES-1))' GES GES xVec(1)]';
            ypatch = [pVec(1:(iGES-1))' yGES 0 0]';
            patch(xpatch, ypatch, GEScolor, 'FaceAlpha', 0.35, ...
                'HandleVisibility', 'off');
        end
        XL(1) = min(min(XL(1), GES-xMargin), xVec(1));
        XL(2) = max(max(XL(2), GES+xMargin), xVec(end));
        xlim(XL);
    end
    xlabel('Chlorophyll \it{a}\rm')
    ylabel('p(Chlorophyll \it{a}\rm)')
end
function [ax, isnew] = checkPlotAx(inpH)
    if(isempty(inpH))
        figH = nFig();
        ax = axes(figH);
        isnew = true;
    else
        ax = inpH;
        isnew = false;
    end
end
%% ---------------
function [ax, C] = getPlotColor(ax)
    % Use axes color order if no color specified
    if(ax.ColorOrderIndex==1)
        ax.ColorOrder = turbo(20);
        newPlot = true;
    end
    C = ax.ColorOrder(ax.ColorOrderIndex,:);
end
%% -----------
function [figH] = nFig()
    targetWidth = 3.25;
    wToH = 2.3;
    pxPerInch = 200;
    figH = figure("Units","pixels",...
        "Position", ...
        [1920 700 pxPerInch*targetWidth pxPerInch*targetWidth/wToH]);
end