function [optExport] = exportOptions(options)
    %EXPORTOPTIONS Set options for exporting figures.
    %   optExport = exportOptions() Get default options for exporting a figure.
    %
    %   optExport = exportOptions(..., 'fname', 'figure_name') Sets file name
    %   (before file extension) of the export figure to 'figure_name'.
    %
    %   optExport = exportOptions(..., 'addTimeStamp', true/false) Adds a 
    %   timestamp to the beginning of the filename.
    %
    %   optExport = exportOptions(..., 'ext', 'png') Sets the file format to
    %   .png. The available formats are 'png', 'eps', 'pdf', 'tiff', 'jpg' and
    %   'gif'. Multile formats can be given as a string array, e.g. ["png",
    %   "gif"], in which case a file for each format is created.
    %
    %   optExport = exportOptions(..., 'dir', pathToFolder) Sets the folder
    %   where the figure will be export to.
    %
    %   optExport = exportOptions(..., 'Resolution', 600) Sets the resolution
    %   (dpi) of the saved figure to 600 (default).
    %
    %   See also exportFig.
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        options.fname {mustBeTextScalar} = 'figure';
        options.addTimeStamp {mustBeNumericOrLogical} = true;
        options.ext {mustBeSubset(options.ext, ...
            {'png', 'eps', 'pdf', 'tiff', 'jpg', 'gif'})} = 'png';
        options.dir {mustBeFolder} = 'figures/';
        options.Resolution (1,1) {mustBePositive} = 600;
    end

    optExport = options;
end

%% --------- validation functions ---
function mustBeSubset(a, b)
    % test that each element in a is in b
    mustBeText(a);
    if iscell(a)
        % multiple elements
        for i = 1:length(a)
            mustBeMember(a{i}, b);
        end
    else
        mustBeMember(a, b);
    end

end