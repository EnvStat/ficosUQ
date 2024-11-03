function [tb] = squishTbVals(tb, varargin)
%SQUISHTBVALS tb = squishTbVals(tb, bounds, valueColumn, indexColumn,
%variable)
%   Detailed explanation goes here

ip = inputParser;

% Boundary, either only lower bound as scalar or 2 element vector
% [lower upper]
tbVars = tb.Properties.VariableNames;
validBoundary = @(x) isnumeric(x) && (isscalar(x) || ...
    length(x)==2 && x(2)>x(1));
ip.addOptional('boundary', 0, validBoundary);

% Variable/column to squish - squish all numeric columns if not given
validColumn = @(x) isempty(x) || all(ismember(x, tbVars));
ip.addOptional('valueColumn', [], validColumn);
% For stacked tables:
% Index variable for stacked variable
ip.addOptional('index', [], validColumn);
% If given, squish values from rows where index == var
validVar = @(x) ischar(x) || ...
    (iscell(x) && all(cellfun(@ischar, x)));
ip.addOptional('var', '', validVar);

parse(ip, varargin{:});

lb = ip.Results.boundary(1);
ub = Inf;
if length(ip.Results.boundary)==2
    ub = ip.Results.boundary(2);
end



if isempty(ip.Results.valueColumn)
    % Squish all numeric columns
    valCol = tbVars(cellfun(@(x) isnumeric(tb.(x)), tbVars));
elseif iscell(ip.Results.valueColumn)
    valCol = ip.Results.valueColumn;
elseif ischar(ip.Results.valueColumn)
    valCol = {ip.Results.valueColumn};
end

squishIdx = true(size(tb,1), 1);
if ~isempty(ip.Results.index) && ~isempty(ip.Results.var)
    squishIdx = ismember(tb.(ip.Results.index), ip.Results.var);
end

for ii = 1:length(valCol)
    col = valCol{ii};
    if lb > -Inf
        tb.(col)(squishIdx) = max(tb.(col)(squishIdx), lb);
    end
    
    if ub < Inf
        tb.(col)(squishIdx) = min(tb.(col)(squishIdx), ub);
    end
end

end

