function [xmode] = gp_mode(gp, xnm, a, ynm, optacq,optimf)
% GP_MODE Locates GP mode.
% [xmode] = gp_mode(gp, xnm, a, ynm, optacq, optimf) Locates the mode of the
% Gaussian Process gp.
%
% Arguments:
% xnm    : training locations
% ynm    : log posterior density at training locations
% a      : vector a = C\ynm, where C is the gp's covariance matrix for xnm 
% optacq : options for the optimization function used for optimisation function
% below.
% optimf : optimisation function used
%
% Copyright (c) 2013 - 2018 Jarno Vanhatalo
% Copyright (c) 2017 - 2024 Karel Kaurila
%

fh_e = @(xp) gp_cov(gp,xp,xnm)*a;
fh_g = @(xp) (grad(xp, fh_e));
fh_eg = @(xp) deal(fh_e(xp), fh_g(xp));
xstart = xnm(ynm==min(ynm),:);
if(size(xstart,1)>1) % remove duplicates
    xstart = xstart(1,:);
end

xmode = optimf(fh_eg, xstart, [], [], [], [], -1*ones(1,size(xnm,2)), ...
    ones(1,size(xnm,2)), [], optacq);
    
end
