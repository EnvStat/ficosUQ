function [optCell] = struct2opt(structIn)
    %STRUCT2OPT Convert a scalar structure array to a cell array of 
    % (name, value) argument pairs passable to a function.
    %  
    % See also NAMEDARGS2CELL
    %
    % Adapted from:
    % https://stackoverflow.com/questions/15013026/how-can-i-unpack-a-matlab-structure-into-function-arguments
    %
    % Copyright (c) 2017-2024 Karel Kaurila
    %
    arguments
        structIn (1,1) struct
    end
    if isMATLABReleaseOlderThan('R2019b')
        % built-in function introduced in R2019b
        % use the stackoverflow solution for older releases
        fnames = fieldnames(structIn);
        vals = struct2cell(structIn);
        optCell = [fnames, vals]';
    else
        % use the built-in function
        optCell = namedargs2cell(structIn);
    end
    optCell = optCell(:);
end

