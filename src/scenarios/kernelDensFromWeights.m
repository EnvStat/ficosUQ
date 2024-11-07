function [P, X] = kernelDensFromWeights(x,w, varargin)
% KERNELDENSFROMWEIGHTS Kernel density estimation with weighted samples.
% [P, X] = kernelDensFromWeights(x, w, ...) Kernel density estimate for samples x
% with weights w. Adapted from 'kernel1' (Copyright (c) 1999 Simo Särkkä).
% 
% Note that here x is a vector, while 'kernel1' returns estimates for each 
% column of a matrix.
%
% P contains smoothed and normalized densities corresponding to coordinates 
% in X. 
%
% See also KERNEL1
%
% Copyright (c) 2017 - 2024 Karel Kaurila
%  

n = length(w);
% Default values in kernel1
% Estimate std using weights
weightedVar = var(x,w);
defS2 = sqrt(weightedVar)/2;

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