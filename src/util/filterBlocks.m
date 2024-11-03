function [simOut] = filterBlocks(simOut, blocks)
    %FILTERBLOCKS Filter simulator output to only given water formations (called
    %blocks internally in FICOS).
    % simOut2 = filterBlocks(simOut, blocks)
    %
    %
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %

    timeVar = simOut.Properties.DimensionNames{1};

    simVars = simOut.Properties.VariableNames;

    blockVar = 'polyID';
    if any(contains(simVars, 'polyID'))
        blockVar = 'polyID';
    elseif any(contains(simVars, 'wf_id'))
        blockVar = 'wf_id';
    elseif any(contains(simVars, 'block'))
        blockVar = 'block';
    end
    simOut = convertvars(simOut, blockVar, 'categorical');
    simOut = groupfilter(simOut, timeVar, @(x) ismember(x, blocks), blockVar);
end

