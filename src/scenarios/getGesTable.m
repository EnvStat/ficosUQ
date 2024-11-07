function [T_GES] = getGesTable()
%GETGESTABLE Criteria for Good Environmental Status.
%
% T_GES = getGesTable() Construct a table with criteria for Good Environmental Status for
% different regions of the Finnish Archipelago sea. Column 'GES value' contains
% threshold for chlorophyll concentration for each region based on European
% Union's Water Framework Directive.
%
% Copyright (c) 2017 - 2024 Karel Kaurila
%
T_GES = table(categorical( ...
    {'Inner Archipelago';'Middle Archipelago'; 'Outer Archipelago'}),...
    [3.0; 2.5; 2.3], ...
    'VariableNames', {'region', 'GES value'});

end

