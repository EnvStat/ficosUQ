function [sOut] = assignFields(sIn, sIn2)
    %ASSIGNFIELDS Assigns values to matching fields in struct sIn from sIn2.
    %   sOut = assignFields(sIn, sIn2)
    %
    %   Used for filling in multiple function arguments from a struct.
    %   
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        sIn (1,1) struct
        sIn2 (1,1) struct
    end
    sOut = sIn;
    commonFields = intersect(fieldnames(sIn), fieldnames(sIn2));
    for fName = commonFields(:)'
        sOut = setfield(sOut, fName{:}, getfield(sIn2, fName{:}));
    end
end

