% GP_MCMC_APPROXIMATION Draw samples from GP emulator using MCMC.
% 
% gp_mcmc_approximation Script for sampling from the GP emulator using MCMC.
%
% Copyright (c) 2013 - 2018 Jarno Vanhatalo
% Copyright (c) 2017 - 2024 Karel Kaurila
%
[xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(...
    x,y,gp,opt,params.transform,lb,ub);

fh_e = @(xp) gp_cov(gp,xp,xnm)*a;
fh_g = @(xp) (grad(xp, fh_e));
fh_eg = @(xp) deal(fh_e(xp), fh_g(xp));

[ynmMin, iMin] = min(ynm);
xstart = xnm(iMin,:);
if(size(xstart,1)>1)
    xstart = xstart(1,:);
end

% Locate GP mode
fprintf('GP mode:')
xmodeGP = gp_mode(gp, xnm, a, ynm, optacq,optimf)
emodeGP = fh_e(xmodeGP)

fprintf('Observed mode:')
xmode = xstart
emode = fh_e(xmode)

% Softer energy function without GP mean for more stable sampling
fh_e = @(x) gp_pred(gp,xnm,ynm,x).*stdynm;

% Using gpstuff implementation 'sls' in gpstuff/mc
gpstuffroot = getGpstuffroot();
mcPath = fullfile(gpstuffroot, 'mc');
oldpath = addpath(mcPath);

try
    [ws,energies,diagn] = sls(fh_e,xstart,params.optsls);
    % Add the mean back in after sampling
    energies = energies + mynm;

    ws_rt = retransform_xy(ws,lb,ub,params.retransform);
catch ME
    disp(ME)
end

path(oldpath);
% Save progress after Slice sampling
save(savpath)