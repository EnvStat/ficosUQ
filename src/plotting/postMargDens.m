function [figParamPost, tlParamPost] = postMargDens(xSamp, params, ...
        w, options)
%POSTMARGDENS Posterior marginal density plots
% For earlier smaller figures,
%   targetWidth = 3.25, pxPerInch = 96, widthToHeight = 2.9
% For current larger, wider figures,
%   targetWidth = 5, pxPerInch = 150, widthToHeight = 4.5
arguments
    xSamp {mustBeNumeric}
    params struct
    w {mustBeNumeric} = [];
    options.sigma {mustBeNumeric} = [];
    options.bins = [];
    options.extra = [];
    options.units (1,1) {mustBeNumericOrLogical} = false;
    options.targetWidth (1,1) {mustBePositive} = 5;
    options.pxPerInch (1,1) {mustBePositive} = 150;
    options.widthToHeight (1,1) {mustBePositive} = 4.5;
    options.ylab {mustBeTextScalar} = 'Posterior density';
    % Compare uniform and weighted density estimates
    options.weightComparison {mustBeNumericOrLogical} = false;
    % Method for estimating density
    options.densMethod {mustBeMember(options.densMethod, ...
        {'kernel'})} = 'kernel';
    % Method for handling importance weights
    options.weightMethod {mustBeMember(options.weightMethod, ...
        {'direct', 'bootstrap'})} = 'direct';
    options.nBootstrap (1,1) {mustBePositive} = 1e4;
end

sigma = options.sigma;
bins = options.bins;
extra = options.extra;

optKernel = {};

optLwd = {'LineWidth',1.2};
if(~isempty(sigma))
    optKernel = [optKernel, {'sigma', sigma}];
end
if(~isempty(bins))
    optKernel = [optKernel, {'bins', bins}];
end
if(~isempty(extra))
    optKernel = [optKernel, {'extra', extra}];
end


% Target width is 3.25 inches
% assume 96 pixels per inch
targetWidth = options.targetWidth;
pxPerInch = options.pxPerInch;
widthToHeight = options.widthToHeight;
figPos = [1920+700, 700, ...
    targetWidth*pxPerInch, targetWidth*pxPerInch/widthToHeight];

targetLabFontSizePt = 6;
labelFontMult = 1.2;
optFontSize = {'FontUnits','normalized',...
    'FontSize', ...
    targetLabFontSizePt*targetWidth/(72*widthToHeight*labelFontMult),...
    'LabelFontSizeMultiplier', labelFontMult};
optLab = {};

d = size(xSamp, 2);
mu_prior = params.mu_prior;
sigma_prior = params.sigma_prior;

figParamPost = figure('Position', figPos);
tlParamPost = tiledlayout(figParamPost, 1, d ,'TileSpacing','tight',...
    'Padding','tight');
% prior marginal density
%priorMargDens = @(x,mu,s)lognpdf(x, mu, s);

thetaNames = params.thetaNames;
plotIS = ~isempty(w);

if strcmp(options.densMethod, 'kernel')
    densEstFcn = @(xs, w) kernelF(xs, w, optKernel{:});
end

if plotIS & strcmp(options.weightMethod, 'bootstrap');
    xBS = datasample(xSamp, options.nBootstrap, 1, 'Weights', w);
end

for i = 1:d
    axParamPost = nexttile(tlParamPost);
    %hHist = histogram(axParamPost, )
    % Kernel density estimate for parameter posterior density
    %[pUnifEst, xUnifEst] = kernelF(xSamp(:,i), optKernel{:});%kernel1(xSamp(:,i));
    [pUnifEst, xUnifEst] = densEstFcn(xSamp(:,i), []);

    hold on;
    maxP = max(pUnifEst);
    [minX, maxX] = bounds(xUnifEst);
    if(plotIS)
        % Weight by importance samples
        if strcmp(options.weightMethod, 'direct')
            [pIS, xIS] = densEstFcn(xSamp(:,i), w);
        elseif strcmp(options.weightMethod, 'bootstrap')
            [pIS, xIS] = densEstFcn(xBS(:,i), []);
        end
        % kernelHist(xSamp(:,i),w);
        % plot(xIS,pIS, '-r', optLwd{:});
        [minXis, maxXis] = bounds(xIS);
        [minX, maxX] = bounds([minX, minXis, maxX, maxXis]);
    end
    xlim(axParamPost, 'tight');
    ylim(axParamPost, 'tickaligned');%'tight');

    set(axParamPost, optFontSize{:});
    % Prior density
    % XL = axParamPost.XLim;
    % YL = axParamPost.YLim;
    
    xt = linspace(minX,maxX, 1000);
    % Normalize prior scale
   
    ytPrior = lognpdf(xt, mu_prior(i), sigma_prior(i));
    ytPrior = 0.95*maxP*ytPrior/max(ytPrior);
    
    % Plotting
    plot(xt, ytPrior,'k-',optLwd{:});
    xPost = xUnifEst;
    pPost = pUnifEst;
    if plotIS && ~options.weightComparison
        xPost = xIS;
        pPost = pIS;
    end
    plot(xPost, pPost, '-b', optLwd{:});
    if(plotIS & options.weightComparison)
        plot(xIS,pIS, '-r', optLwd{:});
    end
    
    xlabel( varNameToYlab(thetaNames{i}, options.units), optLab{:});

    yticks(axParamPost, []);
    if isempty(options.ylab) || i>1
        ylabel(axParamPost,'');
    else
        ylabel(axParamPost, options.ylab);
    end

    box(axParamPost, 'off');
end

end

%% ------------------------------------

function [P, X] = kernelF(x,varargin)
% Adapter for kernel1 and kernelHist that allows setting the same
% parametrization for both unweighted and weighted samples
ip = inputParser;

n = length(x);

defBins = max(50, ceil(sqrt(n)));
defSigma = [];
validW = @(x) isempty(x) || (isvector(x) && isnumeric(x) && all(x>=0));
addOptional(ip, 'w', [], validW);
addParameter(ip, 'sigma', defSigma);
addParameter(ip, 'bins', defBins);
addParameter(ip, 'extra', 0.2);

parse(ip, varargin{:});

if(isempty(ip.Results.w))
    % kernel1 should be the same as uniformly weighted kernelHist
  [P, X] = kernel1(x, ip.Results.sigma, ip.Results.bins, ip.Results.extra);
else
    [P, X] = kernelHist(x, ip.Results.w, ...
        ip.Results.sigma, ip.Results.bins, ip.Results.extra);
end

end

% ---------

function [P, X] = bsFcn(x, w)
    % Kernel estimate through bootstrap sampling
end

% ---------------------

function [P, X] = kernelHist(x,w, varargin)
% Kernel density estimate using importance weights as histogram counts
% Following kernel density estimate in kernel1 by Simo Särkkä
% Using default values for the parameters in kernel1

n = length(w);
% Default values in kernel1
% Estimate std using weights
weightedVar = var(x,w);
defS2 = sqrt(weightedVar)/2;%std(x)/2;

defBins = max(50,ceil(sqrt(n)));
ip = inputParser;

addOptional(ip, 'sigma', [], @(x) isnumeric(x) || isempty(x));
addOptional(ip, 'bins', defBins, @isnumeric);
addOptional(ip, 'extra', 0.2, @(x) isscalar(x)&& isnumeric(x));
parse(ip, varargin{:});

extra = ip.Results.extra;
s = ip.Results.sigma;
if(isempty(s))
    s = defS2;
end
bins = ip.Results.bins;

[xSorted, iSort] = sort(x, 'ascend');
wSorted = w(iSort);
mx = xSorted(end);
mn = xSorted(1);
delta = extra*(mx-mn);
mx = mx + delta/2;
mn = mn - delta/2;
dx = (mx - mn) / bins;

% Here instead of computing sample histogram, we use sample weights
edges = mn:dx:mx;
% bin centers
X = edges(1:end-1)+dx/2;
iBins = discretize(xSorted, edges);
imx = max(iBins);
H = zeros(bins,1);
H(1:imx) = accumarray(iBins, wSorted);

% Smooth with convolution with a Gaussian kernel by following the rest of 
% kernel1 by Simo Särkkä
x = 2*(max(X)-min(X))*(0:(2*size(H,1)-1))'/size(H,1);
G = normpdf(x,mean(x),s);
G = ifftshift(G(:));
P = real(ifft(fft(G,size(G,1)) .* fft(H,size(G,1))));
P = P(1:size(H,1));
P = P / (sum(P) * dx);
P(P<0)=0; % Make sure P >= 0


end