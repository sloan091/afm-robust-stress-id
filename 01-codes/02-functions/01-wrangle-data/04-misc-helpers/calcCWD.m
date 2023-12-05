% =========================================================================
% Name   : calcCWD.m
% Author : Brandon Sloan
% Date   : 2/12/22
%
% DESCRIPTION
% Estimates the cumulative water deficit (CWD) using measured evapotranspiration
% (ET) and precipitation (P) data. The method is taken from Wang-Erlandsson
% et al. (2016).
% 
% INPUTS
%   ET  - Measured evapotranspiration [mm/d]
%   P   - Measured precipitation [mm/d]
%
% OUTPUTS
%   CWD - Cumulative or running water deficit [mm/d]
% =========================================================================
function  CWD = calcCWD(ET,P)

    % Calculate running water deficit [mm]
    WD_t = ET - P;
    CWD = zeros(length(WD_t),1);
    for jj = 1:length(WD_t)-1
        CWD(jj + 1,1) = max(0,WD_t(jj + 1) + CWD(jj));
    end

end