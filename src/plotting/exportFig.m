function [filePaths] = exportFig(fh, optExport)
    %EXPORTFIG Save a figure to file.
    %   exportFig() Exports current figure handle to a file. 
    %   
    %   exportFig(fh) Export figure in figure handle fh to a file.
    % 
    %   Set additional options with 'exportOptions' as follows:
    %   exportOpt = exportOptions(...);
    %   exportArgs = namedargs2cell(exportOpt);
    %   exportFig(fh, exportArgs{:});
    %   
    %   See also 'help exportOptions'.
    %
    %  Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        fh (1,1) = gcf;
        optExport.fname {mustBeTextScalar} = 'figure';
        optExport.addTimeStamp {mustBeNumericOrLogical} = true;
        optExport.ext {mustBeSubset(optExport.ext, ...
            {'png','eps','pdf','tiff','jpg', 'gif'})} = 'png';
        optExport.dir {mustBeFolder} = 'figures/';
        optExport.Resolution (1,1) {mustBePositive} = 600;
    end

    %% Input validation
    if(~isgraphics(fh, 'figure') || ~isa(fh,'matlab.ui.Figure'))
        error('Invalid figure handle.')
    end

    %%

    % Save figure in each of the given extensions
    ext = string(optExport.ext);
    fname = optExport.fname;
    if(optExport.addTimeStamp)
        fname = sprintf('%s_%s',getTimeStamp(2),fname);
    end
    fnames = compose('%s.%s', fname, ext(:));
    filePaths = fullfile(optExport.dir, fnames);
    
    for i = 1:length(filePaths)
        exportgraphics(fh, filePaths{i},...
            'Resolution', optExport.Resolution);
    end


end
%% custom validation functions 
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
