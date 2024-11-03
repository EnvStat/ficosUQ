function [x] = quasiRng(varargin)
%QUASIRNG  Generating quasi-random points. 
%  Similar to pseudo random number generators, except the 'generated' points
%  come from a deterministic low discrepancy sets instead. Used in the Quasi 
%  Monte Carlo literature. 
%
%   x = quasiRng(n, d) Generate n points with dimension d.
% 
%   x = quasiRng(n, d, 'Skip', 100) Generate n with dimension d, stating from 
%   the 101th element of the sequence.  
%
%   Stores parameters persistently, such that the next set of points begins from
%   where the previous set ended. Restarts from beginning if point generation
%   Method or dimension is changed.
%
%   x = quasiRng(0) Return current state of the persistent variables in a struct.
%   The state variables are
%   'nSkip' - the number of points generated (or skipped) since the beginning
%   of the sequence.
%   'd' - the dimension of the generated points
%   'Method' - which quasirandom set or method is used to generate the points. Now
%    only 'sobol' is implemented.
%
%   quasiRng(-1) Reset the state of number generation to the beginning of the 
%   sequence. Retains dimension and method if set previously.
%
%   Copyright (c) Karel Kaurila 2017 - 2024
%
%  See also SOBOLSET
%
persistent nSkip
persistent d
persistent Method

ip = inputParser;

defaultArgs = struct('n', 50, 'd', 5, 'Method', 'sobol', 'Skip', NaN, ...
    'RelativeSkip', false);

ip.addOptional('n', defaultArgs.n, @(x)isscalar(x)&& isnumeric(x));
ip.addOptional('d', defaultArgs.d, @(x)isscalar(x)&& isnumeric(x) && all(x>0));
ip.addParameter('Method', defaultArgs.Method,@(x)ismember(x, {'sobol'}));
validSkip = @(x) isscalar(x) && (isnan(x)) || ...
    ( isnumeric(x) && all(x>=0));
ip.addParameter('Skip', defaultArgs.Skip, validSkip)
validRelSkip = @(x) isscalar(x) && (islogical(x) || isnumeric(x));
ip.addParameter('RelativeSkip', defaultArgs.RelativeSkip, validRelSkip)

%% Parse input
parse(ip,varargin{:});

if ip.Results.n==0
    % Return quasi rng state if called with n=0
    if isempty(nSkip)
        % quasiRng has not been called yet, initialize first with default
        % arguments
        nSkip = 0;
        d = defaultArgs.d;
        Method = defaultArgs.Method;
    end
    x = struct('nSkip', nSkip, 'd', d, 'Method', Method);
    return
elseif ip.Results.n < 0
    % Reset generator state without generating any points
    nSkip = 0;
    if isempty(d)
        d = defaultArgs.d;
    end
    if isempty(Method)
        Method = defaultArgs.Method;
    end
    x = struct('nSkip', nSkip, 'd', d, 'Method', Method);
    return
end

if(isempty(Method))
    % Initialize persistent variables
    d = ip.Results.d;
    Method = ip.Results.Method;
elseif(any(~ismember({'Method','d'},ip.UsingDefaults)) && ...
        (~strcmp(Method,ip.Results.Method) || d~=ip.Results.d ) )
        % Method or dimension changed: new point set
        d = ip.Results.d;
        Method = ip.Results.Method;
end

if(isempty(nSkip))
    if(isnan(ip.Results.Skip))
        nSkip = 0;
    else
        nSkip = ip.Results.Skip;
    end
elseif(~isnan(ip.Results.Skip))
    % 
    if(ip.Results.RelativeSkip)
        nSkip = nSkip + ip.Results.Skip;
    else
        nSkip = ip.Results.Skip;
    end
end

%% Return next batch of numbers 
if(strcmp(Method, 'sobol'))
    P = sobolset(d, "Skip", nSkip);
    n = ip.Results.n;
    x = net(P, n);
    nSkip = nSkip+n;
else
    error('Unknown method %s', Method);
end



end

