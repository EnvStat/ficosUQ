function [xDesign, iDesign, Ddesign] = maximinDesign(xCand, nDesign, ...
        xPrev)
    %MAXIMINDESIGN Experimental design that maximizes minimum distance between
    %design points.
    %   [xDesign, iDesign, Ddesign] = maximinDesign(xCand, nDesign, xPrev)
    %   Constructs a design of nDesign points by picking points from xCand.
    %   Maximizes the minimum distance between all points in the combined design
    %   of [xDesing; xPrev].
    %
    %   Copyright (c) 2017 - 2024 Karel Kaurila
    %
    arguments
        xCand {mustBeNumeric}
        nDesign (1,1) {mustBePositive}
        xPrev {mustBeNumeric} = [];
    end

    % Pairwise (euclidean) distances between candidates
    
    nCand = size(xCand, 1);
    
    
    D = pdist2(xCand, [xCand; xPrev]);
    % Mark distances to itself as Inf
    Dcand = D(1:nCand, 1:nCand);
        Dcand(Dcand==0) = Inf;
    D(1:nCand, 1:nCand) = Dcand;
    % Minimum distances from candidates to others
    D = min(D, [], 2);
    % Order by maximin distance
    [Ddesign, iSorted] = sort(D, 'descend');
  
    Ddesign = Ddesign(1:nDesign);
    iDesign = iSorted(1:nDesign);
    xDesign = xCand(iDesign, :);
    
end

