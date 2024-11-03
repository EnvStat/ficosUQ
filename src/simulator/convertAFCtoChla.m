function [chla] = convertAFCtoChla(t, A, FC)
    %CONVERTAFCTOCHLA Estimate chla from algal biomasses A and FC.
    %   chla = convertAFCtoChla(t, A, FC) Estimates chlorophyll a from algal
    %   biomasses A and FC.
    %
    %   Arguments:
    %     t  : datetime column vector, time of observation / prediction for A,
    %     FC and chla.
    %     A  : biomass for algae other than cyanobacteria
    %     FC : biomass for cyanobacteria
    %
    %   Copyright (c) Karel Kaurila 2017 - 2024
    %
    arguments
        t (:,1) datetime
        A (:,:) {mustBeNonnegative, mustHaveEqRows(A,t)}
        FC (:,:) {mustBeNonnegative, mustBeEqualSize(A,FC)}
    end

    chla = zeros(size(A));
    indSummer = month(t)>=6;
    chla(~indSummer,:) = 5.67857/30*A(~indSummer,:);
    chla(indSummer,:) = 5.67857/20*(A(indSummer,:)+FC(indSummer,:));
end
%% ----------
function mustHaveEqRows(a,b)
    % Test for equal size
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:notEqual';
        msg = 'Row size of first input must equal row size of second input.';
        error(eid,msg)
    end
end
% Custom validation function
function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';
        error(eid,msg)
    end
end