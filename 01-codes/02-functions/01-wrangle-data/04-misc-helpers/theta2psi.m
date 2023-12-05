% =========================================================================
% Name   : theta2psi.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function transforms relative soil water saturation into soil
% water potential using the soil water characteristic equations laid out by
% Clapp and Hornberger (1978). 
%
% INPUTS
%   psi_sat - Soil water potential near saturation (MPa)
%   b       - Brooks-Corey soil water retention exponent [-]
%   s       - Relative soil saturation (soil water content/porosity) [-]
%
% OUTPUTS
%   psi_s    - Soil water potential [MPa]
%
% REFERENCES
%   (1) Clapp, R. B., & Hornberger, G. M. (1978). Empirical equations for 
%   some soil hydraulic properties. Water Resources Research, 14(4), 
%   601–604. https://doi.org/10.1029/WR014i004p00601
%   
% =========================================================================

function psi_s = theta2psi(psi_sat,b,s)

psi_s = psi_sat*s.^(-b);

end




