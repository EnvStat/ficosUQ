function [sigmas,gVariable, gSample] = samplePostPredSigmas(...
    F, y, w, obsModel, options)
%SAMPLEPOSTPREDSIGMAS Sample posterior predictive observation error 
% parameter sigmas. 
%   [sigmas, gVar, gSamp] = samplePostPredSigmas(F,y,w) Posterior predictive 
% sampling for observation error parameter sigmas conditional on observations y 
% and posterior predictions F using (optional) sample weights w.
%
% Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    F (:,:) timetable
    y (:,:) timetable
    w (:,1) {mustBeNumeric} = 1;
    obsModel (1,1) struct = getObsModel();
    options.nSamp (1,1) {mustBePositive} = 1;
    options.burnin (1,1) {mustBePositive} = 2;
end

% Match predictions against observations
tbComp = comparisonTable(F, y, obsModel);

% truncate values to simulator lower bounds
tbComp.prediction = max(tbComp.prediction, tbComp.simLB);
tbComp.observation = max(tbComp.observation, tbComp.simLB);

% group rows according to sigma groupings and samples
[gSigma, gVariable, gSample] = findgroups(tbComp.variable, tbComp.sampleID);


sigmaSampFcn = @(yy, ff, lb) {sampleSigma(yy, ff, lb, ...
    nu0=obsModel.nu0, ...
    nSamp=options.nSamp, burnin=options.burnin)};
sigmas = splitapply(sigmaSampFcn, ...
    tbComp.observation, tbComp.prediction, tbComp.censureLB, gSigma);

end
