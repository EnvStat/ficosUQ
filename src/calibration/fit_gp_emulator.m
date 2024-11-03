function [xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(...
    x,y,gp,opt,transform,lb,ub, options)
    % FIT_GP_EMULATOR Fits GP emulator to function evaluations.
    %
    % [xnm, ynm, mynm, stdynm, gp, K, C, invC, a] = fit_gp_emulator(...
    % x,y,gp,opt,transform,lb,ub)
    % Updates the hyperparameters of the GP emulator to fit current function 
    % evalutions. Also returns the transformed and normalized
    % parameterizations and outputs along with covariance information used for
    % prediction. Based on a recurring code block in the original calibration 
    % script.
    % 
    % Arguments:
    % - x : (n x d) matrix of parametrisations for the function evaluations
    % - y : (n x 1) vector of outputs to the function evaluations corresponding to
    % x
    % - opt : structure containing options for GP hyperparameter optimization.
    % - transform : function handle used for transforming parameters x.
    % - lb, ub : lower and upper bounds for x, respectively.
    %
    % Optional (name, value) arguments:
    % ... = fit_gp_emulator(..., 'yMaxDiff', ymd) Set the treshold value for
    % truncating large y values to ymd, i.e. values y > min(y)+ymd will be 
    % set to ymd.
    %
    % Copyright (c) 2013 - 2018 Jarno Vanhatalo
    % Copyright (c) 2017 - 2024 Karel Kaurila
    arguments
        x {mustBeReal}
        y (:,1) {mustBeReal}
        gp struct
        opt struct
        transform {mustBeA(transform, 'function_handle')}
        lb {mustBeReal}
        ub {mustBeReal}
        options.yMaxDiff (1,1) {mustBeNumeric, mustBePositive} = 10;
    end

    yMaxDiff = options.yMaxDiff;

    if(isfield(gp,'deriv'))
        tmp = x(:,1:(size(x,2)-1));
        xnm(:,1:size(x,2)-1) = 2*(transform(x)-lb')./(ub'-lb')-1;
    else
        xnm = 2*(transform(x)-lb')./(ub'-lb')-1;
    end

    % Truncate large (log) negative log densities y to be at most
    % min(y) + yMaxDiff
    ynm = y;
    ynm(ynm>min(ynm)+yMaxDiff) = min(ynm)+yMaxDiff;

    % standardize y
    mynm = mean(ynm); stdynm = std(ynm);
    if std(ynm)==0
        stdynm = 1;
    end
    ynm = (ynm-mynm)./stdynm;

    gp = gp_optim(gp,xnm,ynm,'opt',opt);

    [K, C] = gp_trcov(gp,xnm);
    % Find the mode of the GP model for the log minus log posterior
    invC = inv(C);
    a = C\ynm;
end
