function [] = writeParamsToSettings(varargin)
    %WRITEPARAMSTOSETTINGS Writes simulator parameters to ficossettings.ini.
    %   writeParamsToSettings(settingsFile, theta, thetaNames) Writes parameters
    %   theta with names thetaNames to given file settingsFile.
    %
    %   Internal function when running the FICOS simulator from Matlab.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %

    ip = inputParser;

    % either an existing file or char vector with path to a new file
    ip.addRequired('settingsFile', @(x)isfile(x)||ischar(x));
    ip.addRequired('theta', @isnumeric);
    ip.addRequired('thetaNames', @iscell);
    % Original ficossettings file used as template
    ip.addOptional('source', 'ficossettings.ini', @isfile);

    parse(ip, varargin{:});

    settingsFile = ip.Results.settingsFile;
    sourceFile = ip.Results.source;
    theta = ip.Results.theta(:)';
    thetaNames = ip.Results.thetaNames(:)';

    if(length(theta)~=length(thetaNames))
        error('Incorrect number of parameters\n. %d parameter names, but %d parameterss.\n',...
            length(thetaNames), length(theta));
    end


    if(~isfile(settingsFile))
        copyfile(sourceFile, settingsFile);
    end

    % Handle special parameters
    if(contains(thetaNames, 'Rmax'))
        % Special parameter defined as
        % RAmax  = Rmax;
        % RFCmax = Rmax;
        iC = contains(thetaNames, 'Rmax');
        thetaNames(iC) = [];
        tmp = theta(find(iC,1,'first'));
        theta(iC) = [];

        theta = [theta tmp tmp];
        thetaNames = [thetaNames {'RAmax', 'RFCmax'}];
    end

    for ii = 1:length(thetaNames)
        paramName = thetaNames{ii};
        cmd = sprintf('sed -i ''/%s/c\\%s %.9g\'' %s',...
            paramName, paramName, theta(ii), settingsFile);
        [status, msg] = system(cmd);
        if(~status)
            error(msg);
        end
    end
end