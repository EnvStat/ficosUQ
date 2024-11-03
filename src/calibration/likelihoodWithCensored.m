function [nll, nlj, sumTb] = likelihoodWithCensored(obs, pred, options)
%LIKELIHOODWITHCENSORED Likelihood function that includes censored observations.
%
%   nll = likelihoodWithCensored(obs, pred) Returns negative log
%   likelihood for prediction pred with observations obs.
%
%   Required arguments:
%   obs - timetable of observations
%   pred - timetable of predictions (simulation output)
%
%   Optional (name,value) arguments:
%   'transform' - Transformation applied to observations and predictions. Can be one
%   of 'log' (default), 'log1p', 'sqrt' or 'none', where 'log' is the natural 
%   logarithm, 'log1p' is the log(1+x) transformation, 'sqrt' is the square root
%   and 'none' is no transformation.
%
%   'logOffset' - Set the offset value 'a' in log(a+x) used for the 'log1p'
%   transformation above. Default 'a'=1.
% 
%   'nu0' - prior degrees of freedom for the multivariate t likelihood. Default
%   nu0 = 4.
%
%   See also: 'load_intens', 'run_sim' and 'readOutput'.
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    obs timetable
    pred timetable
    options.transform {mustBeTextScalar, mustBeMember(options.transform,...
        {'log', 'log1p', 'sqrt', 'none'})} = 'log';
    options.logOffset (1,1) {mustBePositive} = 1;
    options.nu0 (1,1) {mustBePositive} = 4;
    % Developer arguments
    options.logCdfMin (1,1) {mustBeReal} = log(1e-100);
    options.allPoints (1,1) {mustBeNumericOrLogical} = false;
    options.updateTauCensored (1,1) {mustBeNumericOrLogical} = false;
    options.tauOutliers (1,1) {mustBeNumericOrLogical} = true;
    options.censorAFC (1,1) {mustBeNumericOrLogical} = false;
    options.pLowerBound (1,1) {...
        mustBeInRange(options.pLowerBound, 0,1, "exclude-lower")} = 1;
    options.censoredMethod {mustBeTextScalar, mustBeMember(...
        options.censoredMethod, {'normcdf', 'tcdf', 'line'})} = 'tcdf';
    options.censCondSigma {mustBeTextScalar, mustBeMember(...
        options.censCondSigma, {'mean', 'mode', 'scale'})} = 'scale';
end
%% Define fixed parameters

% Column used to determine which variable each row represents
indexVar = 'var';
dataVars = {'DIN1','DIN2','DIP1','DIP2','A','FC','chla'};

% chla is the sum of (scaled) A and FC and is excluded from fitting - 
% otherwise it would be counted twice
predVars = {'DIN1','DIN2','DIP1','DIP2','A','FC'};

joinKeys = {'Aika', 'wf_id', indexVar};

% Censuse lower bounds, before transformations
censureBounds = table( categorical(dataVars)', ...
    [10 7 3 2 0.01 0.5 0.02]', ...
    'VariableNames',{indexVar, 'Lower bound'});
censureBounds = censureBounds(ismember(censureBounds.var, predVars),:);


% Correlation between residues
resCorr = 'none';

%% Input parsing 
trMethod = options.transform;
logOffset = options.logOffset;
% Prior degrees of freedom for sigma_k^2
nu0 = options.nu0;

if(strcmp(trMethod,'log1p'))
    obsTransform = @(y) log(y+logOffset);
    obsRetransform = @(z) max(0,exp(z)-logOffset);
    nlogJacobian = @(z) sum(z);
elseif(strcmp(trMethod,'log'))
    obsTransform = @(y) log(y);
    obsRetransform = @(z) exp(z);
    nlogJacobian = @(z) sum(z);
elseif(strcmp(trMethod,'sqrt'))
    obsTransform = @(y) sqrt(y+1);
    obsRetransform = @(y) y.^2-1;
    nlogJacobian = @(z) sum(log(z));
else
    % No transform
    obsTransform = @(y) y;
    obsRetransform = @(z) z;
    nlogJacobian = @(z) 0;
end

logCdfMin = options.logCdfMin;
% Treat censored observations as point observations at the lower bound instead
% (default false)
allPoints = options.allPoints;
pLowerBound = options.pLowerBound;
updateTauCen = options.updateTauCensored;
censorAFC = options.censorAFC;
tauOutliers = options.tauOutliers;
censoredMethod = options.censoredMethod;
%% Ensure tables have the correct column names

if(ismember('polyID', obs.Properties.VariableNames))
    obs = renamevars(obs, 'polyID', 'wf_id');
    obs.wf_id = categorical(obs.wf_id);
end

if ~strcmp(pred.Properties.DimensionNames{1}, 'Aika')
    pred = cleanSimOutput(pred, 'old');
end
if(ismember('polyID', pred.Properties.VariableNames))
    pred = renamevars(pred, 'polyID', 'wf_id');
    pred.wf_id = categorical(pred.wf_id);
end

%% Match observations to predictions
% Use stacked tables (pivot or melted in other software) and innerjoin

if(~ismember(indexVar, obs.Properties.VariableNames))
    obs = stackVars(obs, predVars);
end
if(~ismember(indexVar, pred.Properties.VariableNames))

    if(ismember(trMethod, {'log','log1p'}))
        % 
        pred = squishSimMins(pred);
    end

    pred = stackVars(pred, predVars);
    pred = renamevars(pred, 'observation', 'prediction');
end

% Take predictions that have matching observations
compTb = innerjoin(obs, pred, 'Keys', joinKeys);
clear('pred', 'obs');


%% Censored values

if(~censorAFC && ~ismember(trMethod, {'log','log1p'}))
    % Truncate A and FC observations when using log transforms, keep as is
    % otherwise
    censureBounds.("Lower bound")(censureBounds.(indexVar)=="A") = 0;
    censureBounds.("Lower bound")(censureBounds.(indexVar)=="FC") = 0;
end

compTb = innerjoin(compTb, censureBounds);
compTb.censored = compTb.observation <= compTb.('Lower bound');

% Assign observation as the lower bound when censored
compTb(compTb.censored, 'observation') = ...
    compTb(compTb.censored, 'Lower bound');
% can be adjusted as proportion of lower bound (default = 1)
compTb{compTb.censored, 'observation'} = pLowerBound * compTb{compTb.censored, 'observation'};



%% Transform values

compTb.prediction = obsTransform(compTb.prediction);
compTb.observation = obsTransform(compTb.observation);

%% Negative log likelihood for each variable separately

% (optional developer output) summary of the fit
sumTb = table();

% Negative log likelihood
nll = 0;
% (Optional developer output) negative log Jacobian
nlj = 0;

sumTbVars = {'var','observation type','nll','nl Jacobian','n',...
    'tauPrior','tau2Prior'...
    'tauPost','tau2Post'};
for i = 1:length(predVars)
    subTb = compTb(compTb.var == predVars{i},:);

    % Decompose observations into non-censored (Obs) and censored (Cen)
    % p(obs, cen) = p(obs)p(cen|obs)
    nObs = sum(~subTb.censored);
    nCen = sum(subTb.censored);
    yObs = subTb.observation(~subTb.censored);
    yCen = subTb.observation(subTb.censored);
    fObs = subTb.prediction(~subTb.censored);
    fCen = subTb.prediction(subTb.censored);
    if(allPoints || ...
            (~censorAFC && ismember(predVars{i},{'A','FC'})) )
        % A and FC are not censored due to lower bounds for measurument
        % sensitivity, these are only truncated to a lower bound when 
        % using log transforms
        nObs = nObs+nCen;
        yObs = [yObs; yCen];
        fObs = [fObs; fCen];
        nCen = 0;
        yCen = [];
        fCen = [];
    end
    %% Negative log likelihood from point (i.e. non-censored) observations 

    % Standardize observations
    tau2_0 = var(yObs);
    tau_0 = sqrt(tau2_0);

    resObs = (yObs-fObs);
    
    zObs = resObs/tau_0;
    
    % Calculate negative lpdf directly, since built-in matlab function
    % applies exp(.) at the end
    % (see 'type mvtpdf' for details)
    nllUnnorm = (nu0+nObs)/2*log(1+sum(zObs.^2)./nu0);
    % log normalization constant, assuming uncorrelated observations with
    % the same scale
    nllDenom = nObs/2*(log(tau2_0)+log(nu0*pi));
    nllGamma = gammaln(nu0/2) - gammaln((nu0+nObs)/2);
    nllObs = nllUnnorm + nllDenom + nllGamma;

    nll = nll + nllObs;
    nljObs = nlogJacobian(yObs);
    nlj = nlj + nljObs;

    % sigma_k^2 | yObs ~ scale-Inv-Chi^2(nu0+nObs,
    % nu0*tau2_0+nuTau2Obs/(nu0+nObs))
    % used for censored observations
    nuTau2Obs = sum( resObs.^2);

    nu = nu0+nObs;
    if(~tauOutliers)
        % (experimental option), not used by default
        % Testing effect of outlier observations 
        iOut = resObs.^2 > 4;
        nuTau2Obs = sum(resObs(~iOut).^2);
        nu = nu0+sum(~iOut);
    end
    tau2 = (nu0*tau2_0+nuTau2Obs)/nu;
    
    tau = sqrt(tau2);

    if nargout>=3
        sumTb = [sumTb; ...
            table(categorical(predVars(i)), categorical({'point'}) ,...
            nllObs, nljObs, nObs, tau_0, tau2_0, tau, tau2, ...
            'VariableNames',sumTbVars)];
    end

    nllCen = 0;
    nljCen = 0;
    
    tau0Cen = tau_0;
    tau2_0Cen = tau2_0;
    tau2Cen = tau2;
    tauCen = tau;

    if(nCen > 0)
        % Censored observations conditional on point observations
        resCen = (yCen-fCen);
        lcdfCen = zeros(size(resCen));
        switch censoredMethod
            case {'normcdf','tcdf'}
                % Point estimate for sigma
                switch options.censCondSigma
                    case 'scale'
                        % posterior scale tau
                        sigmaCond = tau;
                    case 'mean'
                        % posterior mean
                        sigmaCond = sqrt(nu*tau2/(nu-2));
                    case 'mode'
                        % posterior mode
                        sigmaCond = sqrt(nu*tau2/(nu+2));
                end
                zCen = resCen/sigmaCond;
                if strcmp(censoredMethod, 'normcdf')
                    lcdfCen = log(normcdf(zCen));
                else
                    lcdfCen = log(tcdf(zCen, nu));
                end
                lcdfCen = max(lcdfCen, logCdfMin);
                nllCen = - sum(lcdfCen);
            case 'line'
                % product of univariate normals, but marginalize over sigma2
                % with a line integral
                zCen = resCen/tau;
                % Inverse cdf for chi distributed scale
                sInvFcn = @(t)sqrt(chi2inv(t, nu)/nu);
                % compute joint density using log densities, since some
                % probabilities may be small
                expSumLog = @(x) exp(sum(log(x)));
                integrandFcn = @(t) expSumLog(normcdf(sInvFcn(t)*zCen));
                nllCen = -log(integral(integrandFcn,0,1,'ArrayValued',true));
        end
        
        nll = nll+nllCen;

        nljCen = nlogJacobian(yCen);

        if(nargout>= 3 && updateTauCen)
            % (Experimental optional output)
            % Update posterior for tau based on censored observations
            % (default false), i.e. only use point observations to update
            % these
            % Censored observations are investigated conditional on point
            % observations
            tau0Cen = tau;
            tau2_0Cen = tau2;
            
            % Normal approximation to Var[yCen| yObs, f]
            % (Arellano-Valle et al., 2012)
            % in the article censoring threshold is 0, here we have positive
            % threshold c and have to translate
            % c_i = mu_i/sigma in the article becomes c_i = (mu_i-c)/sigma
            % here
            cVec = -resCen;
            c2Vec = cVec.^2;
            varyCen = (c2Vec.*normcdf(cVec)+cVec.*normpdf(cVec)).*...
                (1-normcdf(cVec));
            varyCen = varyCen+...
                (normcdf(cVec)-(cVec-normpdf(cVec)).*normpdf(cVec));
            varyCen = tau2_0Cen*varyCen;

            tau2Cen = (nu*tau2+sum(varyCen))/(nu+nCen);
            tauCen = sqrt(tau2Cen);
        end
    end
    if nargout >= 3
        sumTb = [sumTb; ...
            table(categorical(predVars(i)), categorical({'censored'}) ,...
            nllCen, nljCen, nCen, tau0Cen, tau2_0Cen, tauCen, tau2Cen,...
            'VariableNames',sumTbVars)];
    end
end

end