function [filterArgs] = filterBlocksArg(blocksArg)
%FILTERBLOCKSARG Generate function call argument for filterBlocks.
% Helper function only used internally.
%
% Copyright (c) 2017 - 2024 Karel Kaurila
%


if iscell(blocksArg)
    % multiple sets of blocks
    filterArgs = cellfun(@indBlockArg, blocksArg, 'UniformOutput',false);
else
    filterArgs = indBlockArg(blocksArg);
end


end
%% -----------
function filterArg = indBlockArg(blocksArg)

    blockSets = {'intens', 'Aurajoki'};

    filterArg = categorical();
    if isempty(blocksArg) || ...
            ( isScalarStr(blocksArg) && strcmp(blocksArg,'all'))
        % use all blocks
        filterArg = get_wf_ids('all');
    elseif isScalarStr(blocksArg) && ismember(blocksArg, blockSets)
        blockList = getBlockSet(blocksArg);
        filterArg = blockList(:);
    elseif iscategorical(blocksArg) || isstring(blocksArg)
        % select list individual blocks directly
        if(isstring(blocksArg))
            blocksArg = categorical(blocksArg(:));
        end
        filterArg = blocksArg;
    elseif ischar(blocksArg)
        % select one block with char argument
        filterArg = categorical(string(blocksArg));
    end
end
%% -------
function out = isScalarStr(s)
    out = ischar(s) || (isstring(s) || isscalar(s));
end