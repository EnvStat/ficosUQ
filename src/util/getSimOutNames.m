function [varNames] = getSimOutNames(sublist)
    %GETSIMOUTNAMES List names for simulator outputs as they appear in the 
    % results.hdf5 file.
    %   varNames = getSimOutNames() Lists all variables in simulator outputs.
    %
    %   varNames = getSimOutNames(subsetName) Lists a named subset of the 
    %   simulator outputs.
    %
    %   the named subsets are: 
    %   'calibration' - variables used for calibratin FICOS.
    %   'N2fixFC'     - variables used for estimating N2 fixation.
    %   'scenario'    - variables needed for nutrient load reduction scenarios
    %   'extScenario' - extended list of variables combining 'calibration' and
    %                  'scenario'
    %   
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        sublist (1,:) {mustBeTextScalar} = '';
    end

    switch(sublist)
        case 'calibration'
            varNames = ["cA"
                        "cC"
                        "cDIN_0"
                        "cDIN_1"
                        "cDIP_0"
                        "cDIP_1"];
        case 'N2fixFC'
            varNames = ["cN2fixFC"
                        "height_0"
                        "rN2fixFC"];
        case 'scenario'
            varNames = ["cA"
                        "cC"];
        case 'extScenario'
            varNames = ["cA"
                        "cC"
                        "cDIN_0"
                        "cDIN_1"
                        "cDIP_0"
                        "cDIP_1"
                        "cN2fixFC"
                        "height_0"
                        "rN2fixFC"];
        otherwise
            varNames = ["Cdet_0"
                "Cdet_1"
                "Cv"
                "Ndet_0"
                "Ndet_1"
                "Nv"
                "Pdet_0"
                "Pdet_1"
                "Pfev"
                "Pv"
                "cA"
                "cC"
                "cDIN_0"
                "cDIN_1"
                "cDIP_0"
                "cDIP_1"
                "cN2fixFC"
                "height_0"
                "height_1"
                "rN2fixFC"
                "totN_0"
                "totN_1"
                "totP_0"
                "totP_1"];
    end



end

