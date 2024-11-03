function [tbS] = stackVars(tbU,varargin)
%STACKVARS Converts a table of predictions or observations to the stacked (tall)
%format.
%   tbS = stackVars(tbU) Converts a wide observation table to the tall format.
%
%   The wide observation table has a separate column for each variable 
%   (e.g 'DIN1','DIN2',...), while the tall table has an index column 'var'
%   denoting which variable the observation is for (e.g. 'var'='DIN1') and a
%   column 'observation' with the observed value.
%
%   The stacked (tall) format is more convenient for calculating the likelihood
%   function.
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
%   See also TABLE.STACK TABLE.UNSTACK

ip = inputParser;
ip.addOptional('dataVars', {'DIN1', 'DIN2', 'DIP1', 'DIP2','A','FC','chla'},...
    @iscell);
ip.addParameter('colName', 'observation');
ip.addParameter('idxName', 'var');

parse(ip, varargin{:});

dataVars = ip.Results.dataVars;
colName = ip.Results.colName;
idxName = ip.Results.idxName;

colNames = tbU.Properties.VariableNames;
if ismember(colName, colNames)
    % If already stacked, rename columns instead
    tbS = tbU;
   if ismember('variable', colNames) && ~strcmp(idxName, 'variable')
       tbS = renamevars(tbS, 'variable', idxName);
   end
   return
else
    
    tbS = stack(tbU, dataVars,"IndexVariableName",idxName,...
        'NewDataVariableName',colName);
    % Remove NaN rows
    iNaN = isnan(tbS.(colName));
    tbS(iNaN,:) = [];
end

end

