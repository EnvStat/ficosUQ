function [tbRegion] = getWfWeights(hd_files_path, region_version, options)
%GETWFWEIGHTS Weight water formations by area or volume.
%   tbRegion = getWfWeights() Returns a table with water formation codes divided
%   into three regions (Inner, Middle and Outer Archipelago), weighted by area.
%
%   tbRegion = getWfWeights('weightBy', 'volume') Weight by water formation
%   volume instead. The resulting weights will be very similar to those based on
%   area.
%
%   Copyright (c) 2017 - 2024 Karel Kaurila
%
arguments
    hd_files_path {mustBeFile} = 'data/hd_files.hdf5';
    region_version char = 'Aurajoki';
    options.weightBy {mustBeTextScalar, mustBeMember(options.weightBy,...
        {'area', 'volume'})} = 'area';
end


%%
blocks = h5read(hd_files_path, '/Block_list');
% Use only the studied subset
[wf_inner, wf_mid, wf_outer, wf_regions] = get_wf_ids(region_version);

region = [repmat(categorical("Inner Archipelago"), numel(wf_inner),1); ...
    repmat(categorical("Middle Archipelago"), numel(wf_mid), 1); ...
    repmat(categorical("Outer Archipelago"), numel(wf_outer),1)];
wf_id_region = categorical([wf_inner(:); wf_mid(:); wf_outer(:)]);
tbRegion = table(region, wf_id_region, 'VariableNames', {'region', 'wf_id'});

wf_id = categorical(blocks);
ind_regions = ismember(wf_id,wf_regions);
wf_id = wf_id(ind_regions);
blocks = blocks(ind_regions);

% Read area and volume from hdf5 file
tbDimensions = table(wf_id);
block_area = zeros(length(wf_id),1);
block_vol = zeros(length(wf_id),1);
for i = 1:length(blocks)
    area_path = sprintf('/Blocks/%d/Area/',blocks(i));
    vol_path = sprintf('/Blocks/%d/Volume/',blocks(i));
    block_area(i) = h5read(hd_files_path, area_path, [1 1], [1 1]);
    block_vol(i) = h5read(hd_files_path, vol_path, [1 1], [1 1]);
end

tbDimensions = addvars(tbDimensions, block_area, 'NewVariableNames', 'Area');
tbDimensions = addvars(tbDimensions, block_vol, 'NewVariableNames', 'Volume');

tbRegion = innerjoin(tbRegion, tbDimensions, 'Keys', {'wf_id'});
tbTotals = groupsummary(tbRegion, 'region', 'sum', ["Area", "Volume"]);
tbTotals.GroupCount = [];
tbRegion = innerjoin(tbRegion, tbTotals, 'Keys', {'region'});

switch options.weightBy
    case 'area'
        tbRegion.weight = tbRegion.Area ./ tbRegion.sum_Area;
    case 'volume'
        tbRegion.weight = tbRegion.Volume ./ tbRegion.sum_Volume;
end
tbRegion = tbRegion(:, {'region', 'wf_id', 'weight'});
