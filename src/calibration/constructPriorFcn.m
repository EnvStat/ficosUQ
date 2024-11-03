function [priorFcn, priorMean, priorSigma, limits] = constructPriorFcn(thetaNames)
    %CONSTRUCTPRIORFCN Constructs the parameter prior density function.
    %   priorFcn = constructPriorFcn(thetaNames) Constructs a function handle
    %   for the negative log prior density function for the simulator
    %   parameters. Argument thetaNames is a list of names for the parameters
    %   used. 
    %
    %   [priorFcn, priorMean, priorSigma, limits] = constructPriorFcn(...)
    %   When called with multiple output arguments, this function also returns 
    %   the parameters mu and sigma, as well as the initial lower and upper
    %   bounds for each parameter.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        thetaNames {mustBeText, mustBeMember(thetaNames, ...
            {'Klight', 'LightN2fix', 'LightThres', ...
            'Rmax', 'RAmax', 'RFCmax', 'sedrate'})} = ...
            {'Klight', 'LightN2fix','LightThres', 'RAmax', 'RFCmax'};
    end
    % Define priors with containers.Map, mapping parameter name to vector
    % with prior mean and variance, and initial boundaries
    % Note that prior mean is in log scale,
    % while the limits are in the original scale
    Priors = containers.Map('UniformValues',false);
    Priors('Klight') = [log(10), log(2)/2, 0.5, 40];
    Priors('LightN2fix') = [log(15), log(2)/2, 0.5, 60];
    Priors('LightThres') = [log(10), log(2)/2, 0.5, 40];
    Priors('Rmax') = [-2, 0.4, 0.03, 0.7];
    Priors('RAmax') = [-2, 0.4, 0.03, 0.7];
    Priors('RFCmax') = [-2, 0.4, 0.03, 0.7];
    Priors('sedrate') = [0, log(2)/2, 0.25, 4];

    priorMean = zeros(1, length(thetaNames));
    priorSigma = zeros(1, length(thetaNames));
    limits = zeros(length(thetaNames), 2);
    for i = 1:length(thetaNames)
        prior_vec = Priors(thetaNames{i});
        priorMean(i) = prior_vec(1);
        priorSigma(i) = prior_vec(2);
        limits(i, :) = prior_vec(3:4);
    end
    % negative log density for the log normal prior distribution
    priorFcn = @(x) 0.5*sum( (((log(x)-priorMean)./...
     priorSigma)).^2) + sum(log(x));
end

