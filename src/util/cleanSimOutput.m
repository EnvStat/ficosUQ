function [simOut] = cleanSimOutput(simOut, outFmt)
    %CLEANSIMOUTPUT Convert simulator ouput table into a specific format.
    %   simOut2 = cleanSimOutput(simOut, outFmt) Convert simulator output table into
    %   a specic format required by internal functions.
    %   
    %  Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        simOut timetable
        outFmt {mustBeTextScalar, ...
            mustBeMember(outFmt, {'old','new','mixed'})} = 'new';
    end

    switch outFmt
        case 'new'
            % Newer format
            newTimeVar = 'Time';
            newBlockVar = 'wf_id';
        case 'old'
            % Older format required by some internal functions
            newTimeVar = 'Aika'; % 'Time' in Finnish
            newBlockVar = 'polyID';

        case 'mixed'
            newTimeVar = 'Aika';
            newBlockVar = 'wf_id';
        otherwise
            error('Invalid version %s', outFmt)
    end

    simOut.Properties.DimensionNames{1} = newTimeVar;

    blockVars = {'polyID','wf_id', 'block'};
    varNames = simOut.Properties.VariableNames;
    iBlockVars = ismember(varNames, blockVars);
    if ~any(iBlockVars)
        error('No compatible location column found');
    end
    blockVar = varNames{find(iBlockVars,1,'first')};

    simOut = renamevars(simOut, {blockVar}, {newBlockVar});
    simOut = convertvars(simOut, {newBlockVar}, 'categorical');
end