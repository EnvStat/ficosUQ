function deltaf  = grad(w, func,  varargin)
%GRAD evaluate the gradient with finite differences.
%
%	Description
%       Evaluate the gradient of a function FUNC at a parameter vector X.
%       A central difference  formula with step size 1.0e-6 is used. 
%
%	GRADCHECK(X, GRAD, P1, P2, ...) allows additional arguments to
%	be passed to FUNC and GRAD.
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)
%   Copyright (c) Jarno Vanhatalo (2011)

% This software is distributed under the GNU General Public 
% License (version 2 or later); please refer to the file 
% License.txt, included with the software, for details.


% Reasonable value for step size
epsilon = 1.0e-6;

func = fcnchk(func, length(varargin));

% Treat
nparams = length(w);
deltaf = zeros(1, nparams);
step = zeros(1, nparams);
for i = 1:nparams
  % Move a small way in the ith coordinate of w
  step(i) = 1.0;
  fplus  = feval('linef', epsilon, func, w, step, varargin{:});
  fminus = feval('linef', -epsilon, func, w, step, varargin{:});
  % Use central difference formula for approximation
  deltaf(i) = 0.5*(fplus - fminus)/epsilon;
  step(i) = 0.0;
end

function y = linef(lambda, fn, x, d, varargin)
%LINEF	Calculate function value along a line.
%
%	Description
%	LINEF(LAMBDA, FN, X, D) calculates the value of the function FN at
%	the point X+LAMBDA*D.  Here X is a row vector and LAMBDA is a scalar.
%
%	LINEF(LAMBDA, FN, X, D, P1, P2, ...) allows additional arguments to
%	be passed to FN().   This function is used for convenience in some of
%	the optimization routines.
%
%	See also
%	GRADCHEK, LINEMIN
%

%	Copyright (c) Christopher M Bishop, Ian T Nabney (1996, 1997)

% This software is distributed under the GNU General Public 
% License (version 2 or later); please refer to the file 
% License.txt, included with the software, for details.

% Check function string
fn = fcnchk(fn, length(varargin));

y = feval(fn, x+lambda.*d, varargin{:});
