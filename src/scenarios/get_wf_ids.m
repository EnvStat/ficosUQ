function [wf_inner, wf_mid, wf_outer, wf_regions] = get_wf_ids(regVersion, options)
    %GET_WF_IDS List FICOS' internal codes for water formations.
    %
    % [wf_inner, wf_mid, wf_outer, wf_all] get_wf_ids() Lists FICOS' internal
    % water formation identifiers codes for sub regions of the Finnish Archipelago Sea.
    %
    %  (Optional) input arguments
    %  regVersion : which set of water formations to return, either 'old','new','all'
    %            or one of their aliases. 'old' is based on an earlier dision, 'new'
    %            is a more refined division based on the Aurajoki estuary, and
    %            'all' returns all of the water formations in FICOS without dividing
    %            them into sub regtions.
    %  'hd_file': path to the hd_files.hdf5 file used to collect the inernal FICOS
    %  water formation codes.
    %
    %  Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        regVersion {mustBeTextScalar, mustBeMember(regVersion, ...
            {'old', '1',...
            'new','Aurajoki','2', ...
            'all'})} = 'Aurajoki';
        options.hd_file char = 'data/hd_files.hdf5';
    end
    

    v1Strs = {'old', '1'};
    v2Strs = {'new','Aurajoki','2'};
    allStrs = {'all'};

    switch regVersion
        case v1Strs
            % Original region definition based on Finnish Environment Institute's
            % identifiers for each water formation, where the following encodings 
            % denote which region a water formation belongs to:
            %  'Ls' -> inner archipelago
            %  'Lv' -> middle archipelago
            %  'Lu' -> outer archipelago

            wf_inner = categorical([1900000011
                1900000012
                1900000013
                1900000054
                1900000055
                1900000056
                1900000057
                1900000058
                1900000059
                1900000060
                1900000062
                1900000063
                1900000064
                1900000065
                1900000066
                1900000067
                1900000068
                1900000069
                1900000070
                1900000071
                1900000072
                1900000073
                1900000075
                1900000076
                1900000077
                1900000078
                1900000079
                1900000080
                1900000081
                1900000082
                1900000083
                1900000084
                1900000085
                1900000086]);

            wf_mid = categorical([1900000094
                1900000095
                1900000096
                1900000097
                1900000098
                1900000099
                1900000101
                1900000102
                1900000103
                1900000104
                1900000105
                1900000106]);

            wf_outer = categorical([1900000016
                1900000017
                1900000087
                1900000088
                1900000089
                1900000090
                1900000091
                1900000092
                1900000093]);
        case v2Strs
            % Refined division for inner and middle arhipelago regions
            % based on the Aurajoki estuary

            wf_inner = categorical([1900000056
                1900000057
                1900000058
                1900000059
                1900000060
                1900000061
                1900000062
                1900000063
                1900000064
                1900000065
                1900000066
                1900000067
                1900000068
                1900000069
                1900000072
                1900000079]);

            wf_mid = categorical([1900000096
                1900000097
                1900000098
                1900000099
                1900000100
                1900000101
                1900000102
                1900000103
                ]);

            % Outer archipelago is defined as earlier
            wf_outer = categorical([1900000016
                1900000017
                1900000087
                1900000088
                1900000089
                1900000090
                1900000091
                1900000092
                1900000093]);

        case allStrs
            % Return all of FICOS's water formation codes
            mustBeFile(options.hd_file);
            wf_regions = categorical( ...
                h5read(options.hd_file, '/Block_list'));
            % FICOS has no region division on this scale
            % returning all codes intead
            wf_inner = wf_regions;
            wf_mid = wf_regions;
            wf_outer = wf_regions;
    end

    wf_regions = unique([wf_inner; wf_mid; wf_outer], 'stable');
end