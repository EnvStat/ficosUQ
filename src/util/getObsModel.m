function [obsModel] = getObsModel(options)
    %GETOBSMODEL Helper function for defining the observation model.
    %   obsModel = getObsModel() Return default observation model. 
    % 
    %   obsModel = getObsModel(...) Define observation model by filling in
    %   named arguments:
    %
    %   'name'      : name for the observation model.
    %   'obsVars'   : observed variables used for fitting the model (calibration)
    %   'plotVars'  : variables that are plotted in figures
    %   'censureLB' : censoring lower bound for observations.
    %   'simLB'     : lower bounds for simulator predictions (e.g. refugee
    %                 values for algal biomasses).
    %   'transFcn'  : function handle for trasforming observed variables (e.g.
    %                 @log for log-normal observation model)
    %   'retransFcn': inverse of the 'transFcn' above
    %   'nu0'       : prior degrees of freedom for observation error variance.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.name (1,:) char = 'mvLogT'
        options.obsVars (1,:) cell = {'DIN1','DIN2','DIP1','DIP2','A','FC'}
        options.plotVars (1,:) cell = {'DIN1','DIN2','DIP1','DIP2','chla','A','FC'}
        options.censureLB (1,:) {mustBeNumeric} = [10 10 3 2 0 0]
        options.simLB (1,:) {mustBeNumeric} = [0.7 0.7 0.3 0.3 0.01 0.5]
        options.transFcn (1,1) function_handle = @log
        options.retransFcn (1,1) function_handle = @exp
        options.nu0 (1,1) {mustBeNumeric} = 4
    end

    obsModel = options;

end

