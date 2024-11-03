% GPSTUFF_INSTALL Install gpstuff submodule
%
% Launches gpstuff's own matlab_install script. This should be run from the 
% project root folder (ficosUQ) or preferrably through the main installation script
% ('ficosUQ/install.sh').
%
% Copyright (c) 2017 - 2024 Karel Kaurila
% See also src/submodules/gpstuff/matlab_install.m
%
rootDirName = 'ficosUQ';
currPath = pwd;
[~, tmp] = fileparts(currPath);
if ~strcmp(tmp, rootDirName)
    error('Run this script from %s\n now in: %s', rootDirName, currPath);
end

gpstuffPath = 'src/submodules/gpstuff/'; 
cd(gpstuffPath);
matlab_install;
cd(currPath);


