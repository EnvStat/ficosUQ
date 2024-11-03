function [x, y] = retransform_xy(xnm,lb,ub,retransform,ynm,mynm,stdynm)
% RETRANSFORM_XY Retransforms simulator parameters and corresponding log densities
% to their original scale
%
% [x, y] = retransform_xy(xnm, lb, ub, retransform, ynm, mynm, stdynm)
% Retransform transformed and normalized simulator parameters xnm and standardized 
% log densities ynm into their original scales, x and y, respectively.
% 
% x = retransform_xy(xnm, lb, ub, retransform) Transform only xnm to its
% original scale.
%
% Arguments:
% xnm - matrix of transformed and normalized simulator parameters
% lb, ub - lower and upper bounds for xnm in the transformed scale
% retransform - handle to the inverse function transforming xnm back into x
% ynm - standardized log densities
% mynm - mean of the log densities in original scale
% stdynm - standard deviation of the log densities evaluations in original scale
%
% This function is based on a recurring block of code in the original calibration 
% script. See also fit_gp_emulator.
%
% Copyright (c) 2013 - 2018 Jarno Vanhatalo
% Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    xnm {mustBeReal}
    lb (:,1) {mustBeReal}
    ub (:,1) {mustBeReal}
    retransform {mustBeA(retransform, 'function_handle')} = @(x)exp(x);
    % Optional
    ynm (:,1) {mustBeReal} = NaN(size(xnm,1),1);
    mynm (1,1) {mustBeReal} = NaN(1,1);
    stdynm (1,1) {mustBeReal} = NaN(1,1);
end

   x = retransform((xnm+1)/2.*(ub'-lb')+lb');

   if nargout >=2
       % validate ynm, mynm and stdynm only if needed
       mustBeNonNan(ynm);
       mustBeNonNan(mynm);
       mustBeNonNan(stdynm);
       y = ynm.*stdynm+mynm;
   end
end
