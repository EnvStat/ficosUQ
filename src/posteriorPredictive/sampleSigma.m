function [sigmas, yLatent] = sampleSigma(y, f, lb, prior, options)
    % SAMPLESIGMA Gibbs sampling for parameter sigma.
    %   sigmas = sampleSigma(y, f) Sampling for standard deviation
    %   parameter sigma using Gibbs sampling conditional on observations y and
    %   predictions f. 
    %
    %   [sigmas, yLatent] = sampleSigma(y,f,lb) Gibbs sampling from a truncated 
    %   normal distribution with lower bound lb. Return (optional) matrix
    %   yLatent, where rows for censored observations (y(i) < lb) contain 
    %   samples for the latent observation, and for non-censored observations
    %   yLatent(i,:) = y(i).
    %
    %   Additional arguments:
    %   'nu0'    : prior degrees of freedom for sigma
    %   'kappa0' : prior scale squared for sigma. If kappa0=0, use
    %              scale-invariant prior.
    %   'burnin' : amount of burn-in samples 
    %   
    %   Copyright (c) Karel Kaurila 2017 - 2024
    %
    arguments (Input)
        y (:,1) {mustBeNumeric}
        f (:,1) {mustBeNumeric}
        lb (:,1) {mustBeNumeric} = -Inf;
        prior.nu0 (:,1) {mustBePositive} = 4;
        prior.kappa0 (:,1) {mustBeNonnegative} = 0;
        options.nSamp (1,1) {mustBePositive} = 1;
        options.burnin (1,1) {mustBeNonnegative} = 0;
    end
    arguments (Output)
        sigmas (:,1) {mustBePositive}
        yLatent (:,:) {mustBeNumeric}
    end
    %% Parse optional arguments
    nSamp = options.nSamp;
    burnin = options.burnin;
    nu0 = prior.nu0(1);
    if prior.kappa0(1) == 0
        % scale-invariant prior
        kappa0 = var(y);
    else
        kappa0 = prior.kappa0(1);
    end    
    %% Initialize
    
    nInit = 1 + burnin;
    nTotal = nInit+nSamp;
    yLatent = repmat(y, 1, nTotal);
    sigmas = ones(nTotal,1);
    if isscalar(lb)
        % Use same lower bound for all observations
        lb = lb*ones(size(y));
    end
    nu = nu0(1) + length(y);
    censored = y < lb;
    indCensored = find(censored);
    nCensored = sum(censored);
    kappa = kappa0*ones(nTotal,1);
    
    % Inverse CDF sampling for yLatent | sigma_i, y < lb
    muVec = f(censored);
    yLatFcn = @(sig) muVec + sig*...
        norminv(rand(nCensored,1).*normcdf(lb(censored), muVec, sig));
    %% Gibbs sampling
    for s = 2:nTotal
        res = (yLatent(:,s-1)-f);
        kappa(s) = (nu0*kappa0 + sum(res.^2))/nu;
        sigmas(s) = sqrt(sinvchi2rand(nu, kappa(s)));
        % sample yLatent_s | sigma_s
        yLatent(censored, s) = yLatFcn(sigmas(s));
    end
    % Discard initial value and burnin
    sigmas = sigmas((2+burnin):nTotal);
    if nargout >= 2
        yLatent(:,(2+burnin):nTotal);
    end
end