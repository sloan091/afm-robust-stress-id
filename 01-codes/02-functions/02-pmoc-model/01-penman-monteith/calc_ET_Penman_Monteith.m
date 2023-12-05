% =========================================================================
% Name   : calc_ET_Penman_Monteith.m
% Author : Brandon Sloan
% Date   : 8/20/21
%
% DESCRIPTION
% This script calculates 3 different versions of "Penman-Monteith": 1)
% Priestley-Taylor assuming equilibrium ET, 2) the classic Penman-Monteith,
% and 3) the PM correction by McColl (2020).  This script allows me to pass
% relevant meteorological variables and obtain 3 ET estimate.  I have
% incorporate the canopy scale Medlyn model of Knauer et al. (2018) to
% calculate stomatal conductance.
%
% INPUTS
%   PM_opts - MATLAB structure file containing all required parameters and
%   states listed below in the correct units
%
%   Observations
%       P_a    - Atmospheric pressure [Pa]
%       T_a    - Atmospheric temperature [degrees C]
%       time   - Timestamp of observation as a MATLAB data number
%       dt     - Resolution of observations [hr]
%       G_o    - Minimal canopy conductance [mol air/m^2/s]
%       G_1    - Medlyn slope parameter [kPa^0.5]
%       RH     - Relative humidity in fraction form [-]
%       GPP    - Gross Primary Productivity [micromoles CO2/m^2 GA/s]
%       C_a    - Atmospheric CO2 concentration [micromoles CO2/moles air]
%       z_m    - Measurement height of flux tower [m]
%       h_v    - Vegetation height [m]
%       z_om   - Canopy momentum roughness length [m]
%       U      - Streamwise velocity at measurement height [m/s]
%       u_star - Friction velocity [m/s]
%   P     - Daily precipitation time table [mm]
%   wdw   - Moving window size [days]
%
% OUTPUTS
%   alpha - Mean precipitation depth excluding non-precipitating days [mm]
%   lambda - Rate parameter for waiting time between rainfall events, i.e.,
%   how many rainfall events per day on average [1/day]
%
% =========================================================================
function [LE,ET,VPD_l,In] = calc_ET_Penman_Monteith(Obs,PM_opts)


% PART 1: DEFINE CONSTANTS AND INPUT VARIABLES
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Specific heat at constant pressure [J/kg/K]
c_p = 1004;

% Molecular weight of water [kg/mol]
M_w = 0.018;

% Molecular weight of air [kg/mol]
M_a = 0.029;

% Latent heat of vaporization [J/kg]
L_v = 2.5e6;

% Conversion factor from W/m^2 to mm/day
conv = 86400/L_v;

% Universal gas constant [J/K/mol]
R_g = 8.314;

% Average density of air (kg/m^3)
rho_a = Obs.P_a./(R_g*(Obs.T_a + 273.15))*M_a;


% PART 1: CALCULATE PRIESTLEY TAYLOR ET
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Psychrometric constant [Pa/degrees C]
gma = c_p*Obs.P_a/((M_w/M_a)*L_v);

% Saturated vapor pressure of air [Pa]
e_sat = 611*exp((17.27*Obs.T_a)./(237.3 + Obs.T_a));

% Slope of Clausius-Clapeyron relation [Pa/degrees C]
Delta = (2508.3./(237.3 + Obs.T_a).^2).*e_sat/611*1000;

% Available energy for ET [W/m^2]
Q_ne = Obs.R_n - Obs.G;

% Vapor pressure deficit of air [Pa]
VPD_a = e_sat.*(1 - Obs.RH);

if isequal(PM_opts.CondMethod,1)
    
    % Aerodynamic conductance using log law and Businger-Dyer stability [m/s]
    [~,G_a,~] = calc_g_a(PM_opts.h_v,PM_opts.z_m,PM_opts.z_0m_c,...
        Obs.U,PM_opts.calc_Stability,Obs.T_a,Obs.H,Obs.u_star);
    
elseif isequal(PM_opts.CondMethod,2)
    
    % Empirical aerodynamic conductance as in Knauer (2018) [m/s]
    G_a = 1./(Obs.U./Obs.u_star.^2 + 6.2*Obs.u_star.^0.67);
    
else
    error('Failed to specify valid atmospheric conductance method')
end

% Correct Ga for LAI
if isequal(PM_opts.LAIFlag,1)
    G_a = G_a./Obs.LAI;
else
end

% Conversion from moles air/m^2/s to m/s
c1 = MolFlux2mps(1,Obs.P_a,Obs.T_a);

% Conversion from W/m^2 to moles H2O/m^2/s
c2 = 1/(L_v*M_w);

% Conversion from moles H2O/moles air to kPa
c3 = Obs.P_a/1e3;

% Combine conversion factos
c4  = c2.*c3;

% PART 3: CALCULATE ET
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Penman-Monteith helper fxn [W/m^2]
ET_PM = @(G_c) (Delta.*Q_ne + rho_a.*c_p.*G_a.*VPD_a)./...
    (Delta + gma.*(1 + G_a./(G_c.*c1)));

% Classic Penman-Monteith ET estimate [mm/day]
LE = ET_PM(Obs.G_c);
ET = LE.*conv;

% Leaf level VPD [kPa]
VPD_l = LE.*c4./Obs.G_c;

% Store input variables
ET_i = rho_a.*c_p.*G_a.*VPD_a;
In = table(gma,Delta,Q_ne,ET_i,G_a,c1);

end
