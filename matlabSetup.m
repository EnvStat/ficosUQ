%% Setup script
% Installs Matlab dependencies

% gpstuff
fprintf('Installing gpstuff.\n');
oldpath = addpath('src');
addpath('src/scripts');
gpstuff_install;

path(oldpath);
fprintf('All done.\n');
