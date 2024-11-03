function titleStr = wf_id_to_title(wf_id)
    % WF_ID_TO_TITLE Convert internal water formation code to a title string.
    % titleStr = wf_id_to_title(wf_id)
    %
    % See also WF_ID_TO_STATION
    % 
    % Copyright (c) 2017 - 2024 Karel Kaurila
    %
    switch(categorical(wf_id))
        case '1900000274'
            titleStr = sprintf('Brändö - %s', string(wf_id));
        case '1900000089'
            titleStr = sprintf('Utö - %s', string(wf_id));
        case '1900000099'
            titleStr = sprintf('Seili - %s', string(wf_id));
        otherwise
            titleStr = string(wf_id);
    end
end