function ylabOut = varNameToYlab(varName, units, options)
    % VARNAMETOYLAB Convert internal variable names to axes labels. 
    %
    %  Converts internal variable names for simulator parameters and observed
    %  variables into human readable names used in axes labels. Simulator 
    %  with updated official names will be converted to the newer name.
    %  
    %
    % ylabOut = varNameToYlab(varName) Returns axes label for internal variable.
    %
    % yLabStr = varNameToYlab(varName, units) Add a second row with variable
    % units (units=true) or not (units=false).
    %
    %  Examples: 
    %   >> yl = varNameToYlab('DIN1')
    %      yl = 
    %          {'DIN surface'} {'(\(μg N/L)')
    %
    %
    %   >> varNameToYlab('LightThres', false)
    %      ans = 'ProdThres'
    %
    %   Addtional named arguments:
    %
    %   'compact' (true/false) : use abbreviated names
    %
    %  Copyright (c) 2017 - 2024 Karel Kaurila
    %
    %  See also WF_ID_TO_STATION
    %
    arguments
        varName (1,:) {mustBeTextScalar}
        units (1,1) {mustBeNumericOrLogical} = true;
        % shorter names for article figures
        options.compact {mustBeNumericOrLogical} = false;
    end

    obsNames = {'DIN1', 'DIN2', 'DIP1', 'DIP2', 'A', 'FC', 'chla'};
    unitStr = '';
    varStr = '';
    varStrCompact = '';
    switch string(varName)
        %% Observable names
        case 'DIN1'
            varStr = 'DIN surface';
            varStrCompact = 'DINs';
            unitStr = '(μg N/L)';
        case 'DIN2'
            varStr = 'DIN deep';
            varStrCompact = 'DINd';
            unitStr = '(μg N/L)';
        case 'DIP1'
            varStr = 'DIP surface';
            varStrCompact = 'DIPs';
            unitStr = '(μg P/L)';
        case 'DIP2'
            varStr = 'DIP deep';
            varStrCompact = 'DIPd';
            unitStr = '(μg P/L)';
        case 'A'
            varStr = 'Other algae';
            varStrCompact = 'Algae';
            unitStr = '(μg N/L)';
        case 'FC'
            varStr = 'Cyanobacteria';
            varStrCompact = 'Cyanob.';
            unitStr = '(μg N/L)';
        case 'chla'
            varStr = 'Chlorophyll\it{a}\rm';
            varStrCompact = 'Chl\it{a}\rm';
            unitStr = '(μg\rm/L)';
        %% Parameter names
        case 'Klight'
            varStr = 'Klight';
            varStrCompact = varStr;
            unitStr = '(MJ m^{-2} d^{-1})';
        % Updated names for older internal names    
        case 'LightN2fix'
            varStr = 'N2fixThres';
            varStrCompact = varStr;
            unitStr = '';
        case 'LightThres'
            varStr = 'ProdThres';
            varStrCompact = varStr;
            unitStr = '';
        case 'RAmax'
            varStr = 'RAmax';
            varStrCompact = varStr;
            unitStr = 'd^{-1}';
        case 'RFCmax'
            varStr = 'RFCmax';
            varStrCompact = varStr;
            unitStr = 'd^{-1}';
        otherwise
            varStr = string(varName);
    end
    %
    if isempty(varStrCompact)
        % If compact variable name is not defined, use the full name
        varStrCompact = varStr;
    end
    
    if options.compact
        ylabOut = {varStrCompact};
    else
        ylabOut = {varStr};
    end

    if units
        % add units
        ylabOut = [ylabOut, unitStr];
    end
end