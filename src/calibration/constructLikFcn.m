function [likFcn] = constructLikFcn(tt_data, likName)
    %CONSTRUCTLIKFCN Construct likelihood function.
    %   likFcn = constructLikFcn(tt_data) Constructs a function handle
    %   for the negative log likelihood function. The constructed function 
    %   expects a timetable containing simulator output as its only argument.
    %
    %   Arguments:
    %  
    %   tt_data - (required) timetable containing calibration observations, see
    %   'load_intens'.
    %
    %   likName - (optional) likelihood family used, use tab-completion to see
    %   options. Default 'censoredLog' for log transformed observations with
    %   censoring lower bounds.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        % observation data
        tt_data timetable 
        % name of the likelihood function, see definitions after the arguments
        % block
        likName {mustBeTextScalar, mustBeMember(likName, ...
            {'censoredLog','log', 'censoredLog1p','log1p', ...
            'censoredSqrt','sqrt', ...
            'censored', 'censoredLin', 'censoredNormal'})} = 'censoredLog';
    end
    lik_logNames = {'censoredLog','log'};
    lik_log1pNames = {'censoredLog1p','log1p'};
    lik_sqrtNames = {'censoredSqrt','sqrt'};
    lik_censoredNames = {'censored','censoredLin','censoredNormal'};

    switch likName
        case lik_logNames
            likFcn = @(T)likelihoodWithCensored(tt_data, T, ...
            'transform','log');
        case lik_log1pNames
            likFcn = @(T)likelihoodWithCensored(tt_data, T, ...
            'transform','log1p');
        case lik_sqrtNames
            likFcn = @(T)likelihoodWithCensored(tt_data, T, ...
            'transform','sqrt');
        case lik_censoredNames
            likFcn = @(T)likelihoodWithCensored(tt_data, T, ...
            'transform','none');
        otherwise
            error('Invalid likelihood: %s', likName);
    end
end

