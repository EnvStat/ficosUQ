function [tbBounds] = getCensureBounds(options)
    %GETCENSUREBOUNDS Censoring bounds for observations.
    %   tbBounds = getCensureBounds() Returns a table with censoring lower
    % bounds for each observed variable.
    %
    %   tbBounds = getCensureBounds('ltRefAsCens', true) Treat refugee values
    %   for algal biomasses used in simulator as if they were censoring lower
    %   bounds for these variables.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        % Treat simulator refugee values for algal biomasses as lowerbounds
        options.ltRefAsCens (1,1) {mustBeNumericOrLogical} = false;
    end
    varNames = {'DIN1','DIN2','DIP1','DIP2','A','FC','chla'}';
    cenLB = [10 10 3 2 0 0 0]';

    tbBounds = table(categorical(varNames), cenLB, ...
        'VariableNames',{'variable','censorLB'}, ...
        'RowNames',varNames);
    % refugee values
    tbBounds.refVal = zeros(size(tbBounds,1),1);
    tbBounds{{'A','FC'}, 'refVal'} = [0.01, 0.5]';
    % induced refugee value for chlorohyll:
    refChla = convertAFCtoChla([datetime(2006,1,1); datetime(2006,7,1)], ...
        [0.01; 0.01], [0.5; 0.5]);
    refChla = round(min(refChla),1, 'significant');
    tbBounds{{'chla'}, 'refVal'} = refChla;
    if options.ltRefAsCens
        % Observations below refugee values are considered censored
        tbBounds.censorLB = ...
            max(tbBounds{:, {'censorLB', 'refVal'}}, [], 2);
    end

end