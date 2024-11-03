function [tbOut] = squishSimMins(tbIn)
%ASSIGNSIMMINS Squishes small simulator predictions to their boundaries.
%   tbOut = squishSimMins(tbIn) Sets simulator predictions below their lower
%   boundaries to the lower boundary.
%
%  Copyright (c) Karel Kaurila 2017 - 2024
% 

tbOut = tbIn;

if ismember('variable', tbIn.Properties.VariableNames)
    warning('squishSimMins: input table is already in stacked format.')
    return
end

tbOut.DIN1 = max(tbIn.DIN1, 0.7);
tbOut.DIN2 = max(tbIn.DIN2, 0.7);

tbOut.DIP1 = max(tbIn.DIP1, 0.3);
tbOut.DIP2 = max(tbIn.DIP2, 0.3);

tbOut.A = max(tbIn.A, 0.01);
tbOut.FC = max(tbIn.FC, 0.5);
tbOut.chla = max(tbIn.chla, 0.02);

end

