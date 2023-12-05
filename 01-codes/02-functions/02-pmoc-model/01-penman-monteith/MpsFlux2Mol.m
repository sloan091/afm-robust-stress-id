% =========================================================================
% Name   : MpsFlux2Mol.m
% Author : Brandon Sloan
% Date   : 5/26/21
%
% DESCRIPTION
% This function converts conductance units from mol air/m^2/s to m/s. 
% This conversion factor is mentioned in Oleson et al. (2018) and 
% described in full detail in Bonan (2019).
%
% INPUTS
%   g_mps  - Conductance in velocity units [m/s]  
%   P_atm  - Atmospheric pressure [Pa]
%   T      - Temperature [degrees C]
%
% OUTPUTS
%   g_mol  - Conductance in molar units [mol air/m^2/s]
% 
% REFERENCES
%   (1) Oleson, K. W. et al. (2018). Technical Description of the version 5
%   of the Community Land Model (CLM). 
%
%   (2) Bonan, G. (2019). Climate Change and Terrestrial Ecosystem Modeling.
%   Cambridge University Press. https://doi.org/10.1017/9781107339217
% =========================================================================

function g_mol = MpsFlux2Mol(g_mps,P_atm,T)

% Universal gas constant [J/K/mol]
R_g = 8.314;     

% Calculate conversion factor from Ideal gas law [mol/m^3]
conv = P_atm./(R_g*(T + 273.15));

% Converted conductance [mol air/m^2/s]
g_mol = g_mps.*conv;              

end
