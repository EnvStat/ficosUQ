function [tbChlaSum] = scenario_summaries(tbChlaPred, logSampleWeight, tbGES,...
        options)
    %SCENARIO_SUMMARIES Summary for chlorophyll over growth season for load 
    % reduction scenarios.
    %   [tbChlaSum] = scenario_summaries(tbChlaPred)
    %   Returns a summary table for summer chlorophyll concentration based on
    %   scenario predictions in table tbChlaPred.
    %
    %   [tbChlaSum] = scenario_summaries(tbChlaPred, log_w) Use log
    %   sample weights log_w for each sample. If not given, uniform weights are
    %   used by default.
    %
    %   [tbChlaSum] = scenario_summaries(tbChlaPred, log_w, tbGES)
    %   Define Good Environmental Status thresholds for chlorophyl concentration
    %   with table tbGES.
    %
    %   See also SUMMARISECHLASCENARIOS RUNCATCHMENTSCENARIOS GETGESTABLE
    %
    %   Copyright 2017 - 2024 Karel Kaurila
    %
    arguments
        tbChlaPred {mustBeA(tbChlaPred, {'table', 'timetable'})}
        logSampleWeight {mustBeNumeric} = 0;
        tbGES table = getGesTable()
        options.quantiles (1,1) {mustBeNumericOrLogical} = true;
        options.seasonMonths {mustBeNumeric, mustBeVector, ...
            mustBeInRange(options.seasonMonths, 1, 12)} = (6:9);
        options.tbRegion table = getWfWeights();
    end
    %% Define constants
    
    tbVarNames = @(tb) tb.Properties.VariableNames;

    % grouping variables
    outer_load_groups = {'Internal load', 'Atmospheric load'};
    inner_load_groups = {'Catchment area load', 'Point load'};
    load_groups = [outer_load_groups(:)', inner_load_groups(:)'];
    groupvars_region_year = [{'Year', 'region'}, load_groups(:)', {'SampleID'}];
    groupvars_season = setdiff([groupvars_region_year, 'wf_id'], 'region');

    dataVars = "chla";

    %% Parsing input


    if ismember('sampleID', tbVarNames(tbChlaPred))
        tbChlaPred = renamevars(tbChlaPred, 'sampleID', 'SampleID');
    end

    incQuantiles = options.quantiles;

    % Posterior sample weights
    sampleIDs = sort(unique(tbChlaPred.SampleID));
    w = ones(length(sampleIDs),1);

    if(isscalar(logSampleWeight))
        % Uniform weights
        w = w/length(w);
    else
        % Normalizing log weights
        logSumExp = @(x) max(x)+log(sum(exp(x-max(x))));
        indSamples = ismember(categorical(1:length(logSampleWeight), ...
            1:length(logSampleWeight)),...
            sampleIDs);
        logW_vec = logSampleWeight(indSamples);
        logW_vec = reshape(logW_vec, length(logW_vec),1);
        w = exp(logW_vec - logSumExp(logW_vec));
    end
    T_sample_weights = table(sampleIDs, w, ...
        'VariableNames', {'SampleID', 'SampleWeight'});

    %% Mean over growth season

    if isa(tbChlaPred, 'timetable')
        timeVar = tbChlaPred.Properties.DimensionNames{1};
        inSeason = ismember(month(tbChlaPred.(timeVar)), options.seasonMonths);
        tbChlaPred = tbChlaPred(inSeason,:);
        tbChlaPred.Year = year(tbChlaPred.(timeVar));
        tbChlaPred = groupsummary(tbChlaPred, groupvars_season, 'mean', dataVars);
        meanVars = compose('mean_%s', dataVars);
        tbChlaPred = renamevars(tbChlaPred, meanVars, dataVars);
        tbChlaPred.GroupCount = [];
    elseif ~ismember('Year', tbVarNames(tbChlaPred))
        error(['Provided table must be either a timetable, ' ...
            'or a table containing year means over growth seasons']);
    end

    if ~ismember(tbVarNames(tbChlaPred), 'region')
        % Attach region definition and weights
        tbChlaPred = innerjoin(tbChlaPred, options.tbRegion, 'Keys', 'wf_id');
    end

    %% Group summaries

    % Summaries by region, year and sample
    % Group by year, region, scenario loads and sampleID
    
    % average over wf_ids, weighting by wf_id area
    % area weights previously normalized to sum to one
    G_region_year = groupsummary(tbChlaPred, groupvars_region_year, @(w,x)dot(w,x),...
        {"weight", "chla"});
    % Rename average value to something readable
    G_region_year = renamevars(G_region_year, size(G_region_year,2),'chla');

    % Summaries by region and sample
    % Average over years, uniform weight

    groupvars_region = [{'region'}, load_groups(:)', {'SampleID'}];
    G_year = groupsummary(G_region_year, groupvars_region, "mean", "chla");
    G_year = renamevars(G_year, size(G_year,2), 'chla');

    % Add column of sample Weights
    G_year = innerjoin(G_year, T_sample_weights);

    %% Mean chla by region and scenario
    groupvars_chla = [{'region'}, load_groups(:)'];
    tbChlaSum = groupsummary(G_year, groupvars_chla, ...
        {@(w,x)dot(w,x), @(w,x)var(x,w)}, ...
        {"SampleWeight", "chla"});
    tbChlaSum = renamevars(tbChlaSum, size(tbChlaSum,2)-1+[0, 1], ...
        {'mean(chla)', 'var(chla)'});
    % Include coefficient of variation, since var(chla) tends to be
    % proportional to mean(chla)
    cov_chla = sqrt(tbChlaSum.("var(chla)"))./tbChlaSum.("mean(chla)");
    tbChlaSum = addvars(tbChlaSum, cov_chla, 'NewVariableNames', ...
        'Coefficient of variation (chla)');

    %% Empirical CDF and quantiles

    if(incQuantiles)

        % Empirical CDF: Wihtin each group: Sort by chla, then cumulative sum over
        % weights
        G_cdf = sortrows(G_year, [groupvars_chla, {'chla'}],'ascend');
        groupvars_str = cellfun(@string, groupvars_chla);
        groups_cdf = findgroups(G_cdf(:, groupvars_str));
        cdf_cell = splitapply(@(x){cumsum(x)}, G_cdf.SampleWeight, groups_cdf);
        cdf_vec = vertcat(cdf_cell{:});
        G_cdf = addvars(G_cdf, cdf_vec, 'NewVariableNames', 'CDF');


        % Organize Empirical CDFs as a subtable to be included with other summaries
        % at the end
        % Add auxiliary column sampleRank for indexing as using sampleID for this
        % would shuffle chla and CDF values again
        cdfGroups = [{'region'},load_groups(:)'];
        G_cdf = grouptransform(G_cdf, cdfGroups,...
            @(sID) (1:length(sID))', 'SampleID','ReplaceValues',false);
        G_cdf = renamevars(G_cdf, 'fun_SampleID', 'chlaRank');
        G_cdf.chlaRank = categorical(G_cdf.chlaRank);


        U_CDF = unstack(G_cdf, {'chla', 'CDF'}, 'chlaRank', ...
            'GroupingVariables',cdfGroups);
        U_CDF = mergevars(U_CDF, "chla_x" +digitsPattern, ...
            'NewVariableName', 'chla');
        U_CDF = mergevars(U_CDF, "CDF_x" +digitsPattern, 'NewVariableName','CDF');
        U_CDF = mergevars(U_CDF, {'chla','CDF'}, 'MergeAsTable',true, ...
            'NewVariableName', 'chla CDF');

        % Transpose and pack in a cell array
        fCdfCell = @(T){table(T.('chla')',T.('CDF')',...
            'VariableNames',{'chla','CDF'})};
        CdfCells = rowfun(fCdfCell, U_CDF, 'InputVariables', 'chla CDF', ...
            'OutputVariableNames','chla CDF');

        U_CDF.("chla CDF") = CdfCells{:,1};
        clear('CdfCells');

        % Empirical quantiles
        q_f = @(cdf,x,p) x(find(cdf>=p,1,'first'));
        q_vec = [0.025 0.05 0.15 0.25 0.5 0.75 0.85 0.95 0.975];
        qf_array = cell(size(q_vec));
        for qi = 1:length(q_vec)
            qf_array{qi} = @(cdf, x) q_f(cdf, x, q_vec(qi));
        end

        G_q = groupsummary(G_cdf, groupvars_chla, qf_array, {"CDF", "chla"});

        oldVarNames = compose('fun%d_CDF_chla', 1:length(q_vec));
        newVarNames = compose('q%d',1000*q_vec);
        G_q = renamevars(G_q, oldVarNames, newVarNames);

        % IQR, and 70%, 90% and 95% interval widths
        IQR_vec = G_q.q750-G_q.q250;
        w95_vec = G_q.q975-G_q.q25;
        w90_vec = G_q.q950-G_q.q50;
        w70_vec = G_q.q850-G_q.q150;
        % Quartile coefficient of dispersion (Q3-Q1)/(Q1+Q3)
        qcd_vec = IQR_vec./(G_q.q250+G_q.q750);

        G_IQR = addvars(G_q, IQR_vec, w95_vec, w90_vec, w70_vec, qcd_vec, ...
            'NewVariableNames',{'IQR', 'w95','w90','w70','QCD'});
        G_IQR = renamevars(G_IQR, 'q500', 'median');

        % Join back to other summaries

        rightVarNames = [{'median', 'IQR', 'QCD', 'w95','w90','w70'}, ...
            setdiff(newVarNames, 'q500', 'stable')];
        tbChlaSum = innerjoin(tbChlaSum, G_IQR, 'Keys', groupvars_chla, ...
            'RightVariables', rightVarNames);
        % Include empirical quantiles
        tbChlaSum = innerjoin(tbChlaSum, U_CDF, 'Keys',groupvars_chla, ...
            'RightVariables','chla CDF');
        tbChlaSum = movevars(tbChlaSum, 'chla CDF', 'After', 'median');

    end

    %% Good Environmental Status probabilities
    if ~isempty(tbGES)
        % Probabilities for achieving Good Environmental Status
        % Add column of GES threshold values
        G_year = innerjoin(G_year, tbGES);
        % Is chla below threshold
        within_GES = G_year.chla < G_year.("GES value");
        G_GES = addvars(G_year, within_GES, 'NewVariableNames',"Within GES");
        % Pr(GES)
        groupvars_GES = [{'region'}, ...
            load_groups(:)', {'GES value'}];
        P_GES = groupsummary(G_GES, groupvars_GES, @(w,x)dot(w,x), ...
            {"SampleWeight", "Within GES"});
        P_GES = renamevars(P_GES, size(P_GES,2), "Pr(GES)");
        tbChlaSum = innerjoin(tbChlaSum, P_GES, ...
            'Keys', groupvars_chla, ...
            'RightVariables', {'Pr(GES)','GES value'});
        tbChlaSum = movevars(tbChlaSum, {'GES value', 'Pr(GES)'}, 'After', 'mean(chla)');
    end

end

