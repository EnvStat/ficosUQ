function blockList = getBlockSet(setName)
    % GETBLOCKSET Named list of FICOS's internal water formation codes.
    % blockList = getBlockSet(setName) Returns a list of water formation codes
    % corresponding to the setName.
    %
    % setNames {'Aurajoki','old','1','new', '2'} are used in nutrient load
    % reduction scenarios
    %
    % setNames {'intens', 'intensiiviasemat','calibration'} are all alises to
    % the intensive stations used for calibration.
    %
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %
    if ismember(setName, {'Aurajoki','old','1','new', '2'})
        [~, ~, ~, blockList] = get_wf_ids(setName);
    elseif ismember(setName, {'intens', 'intensiiviasemat','calibration'})
        blockList = categorical([1900000089, 1900000099, 1900000274]);
    else
        error('Unknown block set :%s\n', setName);
    end

end

