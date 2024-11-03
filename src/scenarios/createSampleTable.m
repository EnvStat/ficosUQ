function [tbSamples] = createSampleTable(samples, options)
    %CREATESAMPLETABLE Create table for posterior samples and their log
    %densities.
    %   tbSamples = createSampleTable(samples) Creates a table for posterior
    %   samples. Variable 'samples' is a matrix with samples on rows and 
    %   parameters on columns.
    %
    %   tbSamples = createSampleTable(..., 'energies', energies) Define proposal
    %   energies (log proposal densities) for each sample.
    %
    %   tbSamples = createSampleTable(..., 'logPost', logPost) Define log
    %   posterior densities for each sample.
    %
    %   tbSamples = createSampleTable(..., 'logW', logW) Define log sample
    %   weights directly instead of computing them from the proposal and 
    %   posterior densities.
    %
    %   tbSamples = createSampleTable(..., 'sampleIDs', sampleIDs) Set custom
    %   indetifiers for each sample. By default these will be 1:size(samples,1).
    %
    %   createSampleTable(..., 'savePath', path) Save the resulting table under
    %   given path.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        samples {mustBeReal}
        options.varNames {mustBeText} = '';
        options.energies (:,1) {mustBeReal} = NaN(size(samples,1),1);
        options.logPost  (:,1) {mustBeReal} = NaN(size(samples,1),1);
        options.logW (:,1) {mustBeReal} = NaN(size(samples,1),1);
        options.sampleIDs (:,1) = categorical(1:size(samples,1))';
        options.savePath {mustBeTextScalar} = '';
    end

    if all(isnan(options.logW))
    % Modelling is done with negative log densities, while this function uses
    % log densities. Majority or at least some of them should be negative in
    % these settings, thus flipping sign if all positive.

    % log proposal density
    if all(~isnan(options.energies)) && all(options.energies > 0)
        logQ = -options.energies;
    else
        logQ = options.energies;
    end

    % log posterior density
    if all(~isnan(options.logPost)) && all(options.logPost > 0)
        logP = -options.logPost;
    else
        logP = options.logPost;
    end

    % log importance weights
    if all(~isnan(logQ)) && all(~isnan(logP))
        logW = logP - logQ;
    else
        logW = NaN(size(logQ));
    end

    tbSamples = table(options.sampleIDs, samples, logQ, logP, logW, ...
        'VariableNames',{'sampleID','sample','logQ','logP','logW'});

    else
        % Define log weights directly
        tbSamples = table(options.sampleIDs, samples, options.logW, ...
            'VariableNames',{'sampleID', 'sample', 'logW'});
    end

    if ~isempty(options.savePath) && isfolder(fileparts(options.savePath))
        if isfolder(options.savePath)
            outFile = fullfile(options.savePath, 'sampleTable.mat');
        else
            [fParent,~,fext] = fileparts(options.savePath);
            if strcmp(fext, '.mat')
                outFile = options.savePath;
            else
                outFile = fullfile(fParent, 'sampleTable.mat');
            end
        end
        fprintf('Saving sample table to %s\n', outFile);
        save(outFile, "tbSamples");
    end

end

