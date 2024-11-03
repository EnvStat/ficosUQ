function [tt] = load_intens(incBrando, colNameStyle, SeiliMode, options)
% LOAD_INTENS Load calibration data from intensive stations.
%    tt = load_intens() Load calibration data from intensive stations.
%    
%    Optional arguments used internally by various Matlab scripts.
% 
%    Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    % include data from the Brando station
    incBrando (1,1) {mustBeNumericOrLogical} = false;
    % whether to use old or new table column naming style
    colNameStyle {mustBeTextScalar, ...
        mustBeMember(colNameStyle, {'old', 'new', 'FICOS'})} = 'old';
    % which data to include from the Seili station
    % 'all' - include everything
    % 'summer' - only include data from 1. Jun - 31. Oct each year
    % 'summer-algae' - only include algae data ('A', 'FC', and 'chla') 
    %    from 1. Jun - 31. Oct each year
    SeiliMode {mustBeTextScalar, ...
        mustBeMember(SeiliMode, {'all', 'summer', 'summer-algae'})} = 'all';
    options.dataPath {mustBeFile} = ...
        fullfile('data/intensiveStationsData.xlsx');
    % print information to console
    options.verbose {mustBeNumericOrLogical} = false;
end

%% Data loading

% Load sheets separately and combine

tt = table();
if(incBrando)
    % only load Brando if specified
    tt = [tt; readtable(options.dataPath,'sheet',2)];
end
% load Seili, then Utö
tt_seili = readtable(options.dataPath,'sheet',3);
switch(SeiliMode)
    case 'all'
        if options.verbose
            disp('Using all Seili data');
        end
    case {'summer', 'summer-algae'}
        ind_summer = month(tt_seili.Date)>=6 & month(tt_seili.Date) <=10;
        tt_seili = tt_seili(ind_summer,:);
        if strcmp(SeiliMode, 'summer')
            if options.verbose
                disp('For Seili, only using data from calendar days 1.6. - 31.10.')
            end
        elseif strcmp(SeiliMode, 'summer-algae')
            if options.verbose
                disp('For Seili, only using algae data from calendar days 1.6. - 31.10.')
            end
            tt_seili.cDIN_0(:) = NaN;
            tt_seili.cDIN_1(:) = NaN;
            tt_seili.cDIP_0(:) = NaN;
            tt_seili.cDIP_1(:) = NaN;
        end
    otherwise
        error('Invalid Seili option')
end
% Finally load data from Utö
tt = [tt; tt_seili;
    readtable(options.dataPath,'sheet',4)];

% Convert to timetable
tt = table2timetable(tt);


%% Change column names

% Format table column names in given style
switch colNameStyle
    case 'old'
        % Older variable naming scheme for backwards compatibility
        newNames = {'Block','cDIN_0','cDIN_1','cDIP_0','cDIP_1','cA','cC','Chla'};
        oldNames = {'polyID','DIN1','DIN2','DIP1','DIP2','A','FC','chla'};
        if(isMATLABReleaseOlderThan('R2020a'))
            % Renamevars is only available from R2020a (= v.9.8) onwards
            tt = renamevars(tt, newNames, oldNames);
        else
            % Work around for old Matlab
            for i = 1:length(newNames)
                tt.Properties.VariableNames{newNames{i}} = oldNames{i};
            end
        end
        tt.Properties.DimensionNames{1} = 'Aika'; 
    case {'new', 'FICOS'}
        % Match variable names in FICOS
        tt = convertvars(tt,'Block','categorical');
end

end