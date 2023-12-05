% =========================================================================
% Name   : calc_Penman_Monteith.m
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
%   PM_inp - MATLAB structure file containing all required parameters and
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
function [PET_pt,PET_pt_d,PET_pm,PET_pm_d,PET_mc,PET_mc_d] = ...
    calc_Penman_Monteith(PM_inp)

% PART 1: DEFINE PHYSICAL CONSTANTS
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

% Specific gas constant for water vapor [J/k/kg]
R_v = 461;    

% Average density of air (kg/m^3)
rho_a = 1.2;        
 

% PART 1: CALCULATE PRIESTLEY TAYLOR ET
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Psychrometric constant [Pa/degrees C]
gma = c_p*PM_inp.P_a/((M_w/M_a)*L_v);

% Saturated vapor pressure of air [Pa]
e_sat = 611*exp((17.27*PM_inp.T_a)./(237.3 + PM_inp.T_a));

% Slope of Clausius-Clapeyron relation [Pa/degrees C]
%Delta = 4098*e_sat./(237.3 + PM_inp.T_a).^2;
Delta = (2508.3./(237.3 + PM_inp.T_a).^2).*e_sat/611*1000;

% Available energy for ET [W/m^2]
Q_ne = PM_inp.R_n - PM_inp.G;

% Priestely-Taylor ET estimate [PET_pt; mm/day]
PET_pt = (PM_inp.alpha*Delta./(Delta + gma).*Q_ne)*conv;

% Aggregate PET_pt to daily PET in depth units [mm/day]
PET_pt_d = aggData(PM_inp.time,PET_pt,PM_inp.dt,PM_inp.AggMethod);


% PART 3: CALCULATE CLASSIC PENMAN-MONTEITH ET
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Vapor pressure deficit [kPa]
VPD = e_sat.*(1 - PM_inp.RH)/1000;

% Canopy scale conductance [mol air/m^2 GA/s]
if isequal(PM_inp.G_c_prescribe,1)
    
    G_c = PM_inp.G_c;
    
else
    
    % Canopy scale conductance [mol air/m^2 GA/s]
    G_c = PM_inp.G_o + 1.6*(1 + PM_inp.G_1./sqrt(VPD)).*...
        PM_inp.GPP./PM_inp.C_a;
    
end

if isequal(PM_inp.CondMethod,1)
    
    % Aerodynamic conductance using log law and Businger-Dyer stability [m/s]
    [~,G_a,~] = calc_g_a(PM_inp.h_v,PM_inp.z_m,PM_inp.z_0m_c,...
        PM_inp.U,PM_inp.calc_Stability,PM_inp.T_a,PM_inp.H,PM_inp.u_star);
    
elseif isequal(PM_inp.CondMethod,2)
    
    % Empirical aerodynamic conductance as in Knauer (2018) [m/s]
    G_a = 1./(PM_inp.U./PM_inp.u_star.^2 + 6.2*PM_inp.u_star.^0.67);
    
else
    error('Failed to specify valid atmospheric conductance method')
end

% Convert canopy conductance from molar to volume units [m/s]
G_c_p = MolFlux2mps(G_c,PM_inp.P_a,PM_inp.T_a);

% Aerodynamic and Canopy conductance in series [m/s]
if isinf(G_c)
    G_w = G_a;
elseif isinf(G_a)
    G_w = G_s;
else
    G_w = G_a.*G_c_p./(G_a + G_c_p);
end

% Classic Penman-Monteith ET estimate [PET_pt; mm/day]
PET_pm = (Delta.*Q_ne + rho_a*c_p*G_a.*VPD*1000)./...
    (Delta + gma.*(1 + G_a./G_c_p))*conv;

% Aggregate PET_pm to daily PET in depth units [mm/day]
PET_pm_d = aggData(PM_inp.time,PET_pm,PM_inp.dt,PM_inp.AggMethod);

% PART 3: CALCULATE MCCOLL (2020) PENMAN-MONTEITH ET
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Convert vapor pressure to specific humidity [kg H20/kg air]
q_sat =  e_sat.*((M_w/M_a)./PM_inp.P_a);
q_a = q_sat.*PM_inp.RH;

% Messy Lambert W portion
W_o = lambertw(((L_v^2*q_sat)./(c_p*R_v*(PM_inp.T_a + 273.15).^2).*1./(G_a./G_c_p + 1)).*...
    exp(L_v./(R_v*(PM_inp.T_a + 273.15).^2).*(Q_ne + G_w.*rho_a*L_v.*q_a)./(rho_a*c_p*G_a)));

% McColl (2020) revised Penman-Monteith ET estimate [PET_pt; mm/day]
PET_mc = (rho_a*c_p*G_a.*W_o./(L_v./(R_v*(PM_inp.T_a + 273.15).^2)) -...
    rho_a*L_v.*G_w.*q_a)*conv; 

% Correct numerical issues for the calm limit (G_a goes to 0)
calm = isinf(PET_mc) & Q_ne > 0;
PET_mc(calm) = Q_ne(calm)*conv;
    
% Aggregate PET_mc to daily PET in depth units [mm/day]
PET_mc_d = aggData(PM_inp.time,PET_mc,PM_inp.dt,PM_inp.AggMethod);

end
