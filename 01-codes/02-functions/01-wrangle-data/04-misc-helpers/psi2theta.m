% =========================================================================
% Name   : psi2theta.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function transforms soil water potential into soil water content
% using the soil water characteristic equations laid out by
% Clapp and Hornberger (1978). 
%
% INPUTS
%   theta_sat - Saturated soil water content or porosity [m^3 water/m^3 soil]
%   psi_s     - Soil water potential [MPa]
%   psi_sat   - Soil water potential near saturation (MPa)
%   b         - Brooks-Corey soil water retention exponent [-]
%
% OUTPUTS
%   theta_s - Soil water content [m^3 water/m^3 soil]
%
% REFERENCES
%   (1) Clapp, R. B., & Hornberger, G. M. (1978). Empirical equations for 
%   some soil hydraulic properties. Water Resources Research, 14(4), 
%   601–604. https://doi.org/10.1029/WR014i004p00601
%   
% =========================================================================

function theta_s = psi2theta(theta_sat,psi_s,psi_sat,b)

theta_s = theta_sat.*(psi_s./psi_sat).^(-1./b);

end




