function H = hessian(w0, fh_e, fh_g, varargin)
% HESSIAN compute the hessian matrix using finite difference method

% Copyright (c) 2013 Jarno Vanhatalo

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

m = length(w0);
e0 = feval(fh_e,w0,varargin{:});
delta = 1e-6;
H = -1*ones(m,m);

% Compute first using gradients
% If Hessian is singular try computing with
% larger step-size
if ~isempty(fh_g)
    while any(eig(H)<0) && delta < 1e-2
        for i = 1:m
            for j = i:m
                w1 = w0; w2 = w0;
                w1(j) = w1(j) + delta;
                w2(j) = w2(j) - delta;
                
                g1 = feval(fh_g,w1,varargin{:});
                g2 = feval(fh_g,w2,varargin{:});
                
                H(i,j) = (g1(i)-g2(i))./(2.*delta);
                H(j,i) = H(i,j);
            end
        end
        delta = delta + 1e-3;
   end
end
% If the Hessian is still singular or the delta is too large
% try to compute with finite differences for energies.
if any(eig(H)<0) || delta > 1e-2 || isempty(fh_g)
    delta = 1e-6;
    for i=1:m
        w1 = w0; w4 = w0;
        w1(i) = [w1(i)+2*delta];
        w4(i) = [w4(i)-2*delta];
        
        e1 = feval(fh_e,w1,varargin{:});
        e4 = feval(fh_e,w4,varargin{:});
        
        H(i,i) = (e1 - 2*e0 + e4)./(4.*delta.^2);
        for j = i+1:m
            w1 = w0; w2 = w0; w3 = w0; w4 = w0;
            w1([i j]) = [w1(i)+delta w1(j)+delta];
            w2([i j]) = [w2(i)-delta w2(j)+delta];
            w3([i j]) = [w3(i)+delta w3(j)-delta];
            w4([i j]) = [w4(i)-delta w4(j)-delta];
            
            e1 = feval(fh_e,w1,varargin{:});
            e2 = feval(fh_e,w2,varargin{:});
            e3 = feval(fh_e,w3,varargin{:});
            e4 = feval(fh_e,w4,varargin{:});
            
            H(i,j) = (e1 - e2 - e3 + e4)./(4.*delta.^2);
            H(j,i) = H(i,j);
        end
    end
end