function [timeStampChr] = getTimeStamp(stampType)
    %GETTIMESTAMP Represent current time as a character vector.
    % tmStmp = getTimeStamp() Returns a char array with a timestamp of the
    % current time in format 'dd-MMM-uuuu-HH:mm:ss'.
    %
    % tmStmp = getTimeStamp(2) Returns a timestamp in a more compact format,
    % 'uuuuMMdd'.
    %
    % See also Matlab's built-in datestr and datetime.
    %
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        % which type/format of timestamp:
        %  1 (default) 'dd-MMM-uuuu-HH:mm:ss'
        %  2 (compact) 'uuuuMMdd'
        stampType (1,1) {mustBeInteger} = 1;
    end
    stampFormat = '';
    switch(stampType)
        case 1
            stampFormat = 'dd-MMM-uuuu-HH:mm:ss';
        case 2
            stampFormat = 'uuuuMMdd';
        otherwise
            stampFormat = 'dd-MMM-uuuu-HH:mm:ss';
    end
    timeStampChr = ...
        sprintf('%s', datetime('now', 'Format', stampFormat));
end

