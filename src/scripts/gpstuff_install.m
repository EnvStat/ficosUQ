% GPSTUFF_INSTALL Install gpstuff submodule
%
% Launches gpstuff's own matlab_install script. This should be run from the 
% project root folder or preferrably through the main installation script
% 'install.sh'.
%
% Copyright (c) 2017 - 2024 Karel Kaurila
% See also src/submodules/gpstuff/matlab_install.m
%
currPath = pwd;

gpstuffPath = 'src/submodules/gpstuff/'; 
installerFile = fullfile(gpstuffPath,'matlab_install.m');
if ~isfolder(gpstuffPath)
    error(['gpstuff not found. ;ake sure this script is run from the ' ...
        'project root folder, now in:\n%s'], currPath);
end

if ~isfile(installerFile)
    error("%s not found. Please use the main installation script " + ...
        "'install.sh' for installing dependencies.", installerFile)
end
cd(gpstuffPath);
matlab_install;
cd(currPath);


