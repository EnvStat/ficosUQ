function [logW, W] = postPredImpWeights(tbPred, tbSamples, params)
    %POSTPREDIMPWEIGHTS Importance weights for posterior predictions.
    %   logW = postPredImpWeights() Compute log importance weights for posterior
    %   predictive simulations.
    %
    %   Copyright (c) 2017-2024 Karel Kaurila
    %
    arguments
        tbPred timetable
        tbSamples table
        params struct
    end
    uqSamp = unique(tbPred.sampleID);
    nUqSamp = length(uqSamp);
    tbSamples = tbSamples(ismember(tbSamples.sampleID, uqSamp), :);
    logP = NaN(nUqSamp,1);
    for i = 1:nUqSamp
        idSamp = uqSamp(i);
        xSamp = tbSamples.sample(tbSamples.sampleID==idSamp,:);
        fSamp = tbPred(tbPred.sampleID==idSamp,:);
        fSamp = cleanSimOutput(fSamp, 'old');
        logP(i) = -1*params.post_f(xSamp, fSamp);
    end
    tmp = table(uqSamp(:), logP, 'VariableNames',{'sampleID', 'logP'});
    tbSamples.logP = [];
    tbSamples.logW = [];
    tbSamples = innerjoin(tbSamples, tmp, 'Keys', {'sampleID'});
    logW = tbSamples.logP - tbSamples.logQ;
    W = exp(logW - logSumExp(logW));
end

