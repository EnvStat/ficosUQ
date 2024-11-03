% SET_GP_EMU Initialize Gaussian Process emulator.
% 
% Initializes covariance function for the GP emulator, including priors for
% hyper parameters and settings for optimizing these parameters during
% calibration.
%
% Copyright (c) 2013 - 2018 Jarno Vanhatalo
% Copyright (c) 2017 - 2024 Karel Kaurila
%

coeff_sigma_prior = prior_t;
cfc = gpcf_constant('constSigma2',3,'constSigma2_prior', coeff_sigma_prior);
cfse = gpcf_matern52('lengthScale',1, 'magnSigma2_prior', prior_gaussian,...
    'lengthScale_prior',prior_t);
cfse2 = gpcf_matern52('lengthScale',1*ones(1,params.d),...
    'lengthScale_prior',prior_invt,...
    'magnSigma2_prior', prior_t);
cfl = gpcf_linear('coeffSigma2', 3, 'coeffSigma2_prior', coeff_sigma_prior);
cfl2 = gpcf_squared('coeffSigma2', 3, ...
    'coeffSigma2_prior', coeff_sigma_prior,...
    'interactions', 'off');
lik = lik_gaussian('sigma2', 1e-2, 'sigma2_prior', prior_fixed);

gp = gp_set('cf', {cfc, cfl, cfl2, cfse, cfse2}, 'lik', lik);

%% Options for the optimizer of the acquisition function
optimf = @fmincon;
optdefault=struct('GradObj','on','LargeScale','on',...
    'Algorithm','interior-point','TolFun',1e-10,'TolX',1e-4,...
    'Display', 'off');
optacq=optimset(optdefault);

% options for the GP optimizer
opt=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
