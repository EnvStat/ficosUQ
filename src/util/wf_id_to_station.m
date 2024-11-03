function stationStr = wf_id_to_station(wf_id)
    %WF_ID_TO_STATION Convert internal water formation code wf_id to a more 
    %  readable name station name.
    %
    % Also works with vectorized input.
    %
    % stationName = wf_id_to_station(wf_id)
    %
    % Example:
    % >> wf_id_to_station("1900000089")
    % ans = 
    %       'Utö'
    %
    % Copyright (c) 2017-2024 Karel Kaurila
    %
    % See also VARNAMETOYLAB
    %
    if isscalar(wf_id)
        stationStr = wfId2StationScalar(wf_id);
    else
        stationStr = arrayfun(@(x)wfId2StationScalar(x), wf_id, ...
            'UniformOutput',false);
    end
end
%% ----------
function stationStr = wfId2StationScalar(wf_id)
    % Scalar version of the function
    wf_id = categorical(wf_id);
    switch(wf_id)
        case '1900000274'
            stationStr = 'Brändö';
        case '1900000089'
            stationStr = 'Utö';
        case '1900000099'
            stationStr = 'Seili';
        otherwise
            stationStr = char(wf_id);
    end
end