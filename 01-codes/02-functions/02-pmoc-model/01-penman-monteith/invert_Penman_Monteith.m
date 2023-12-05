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
function [alpha_pt,alpha_pt_d,G_c_pm,G_c_pm_d,VPD_l_pm,G_c_mc,G_c_mc_d,VPD_l_mc,T_s,Omega] = ...
    invert_Penman_Monteith(PM_inp)

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

% Universal gas constant [J/K/mol]
R_g = 8.314; 

% Average density of air (kg/m^3)
rho_a = PM_inp.P_a./(R_g*(PM_inp.T_a + 273.15))*M_a;
 

% PART 2: INVERT PRIESTLEY TAYLOR TO FIND ALPHA
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

% Invert Priestley-Taylor alpha parameter from LE
alpha_pt = PM_inp.LE./(Delta./(Delta + gma).*Q_ne);

% Aggregate PET_pt to daily PET in depth units [mm/day]
alpha_pt_d = aggData(PM_inp.time,alpha_pt,PM_inp.dt,'mean');


% PART 3: INVERT CLASSIC PENMAN-MONTEITH TO FIND CANOPY CONDUCTANCE, G_c
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

% Vapor pressure deficit [Pa]
VPD = e_sat.*(1 - PM_inp.RH);

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

% Invert classic Penman-Monteith to find canopy/surface conductance from LE
% observations [m/s]
G_c_pm_v = (gma.*PM_inp.LE.*G_a)./(Delta.*Q_ne + ...
    rho_a.*c_p.*G_a.*VPD - PM_inp.LE.*(Delta + gma));

% Convert from volume [m/s] to molar units [moles air/m^2 GA/s]
G_c_pm = MpsFlux2Mol(G_c_pm_v,PM_inp.P_a,PM_inp.T_a);

% Aggregate G_c_pm to daily PET in depth units [mm/day]
G_c_pm_d = aggData(PM_inp.time,G_c_pm,PM_inp.dt,'mean');

% Calculate an estimate for VPD at the leaf [Pa]
VPD_l_pm = PM_inp.LE.*gma./(rho_a.*c_p.*G_c_pm_v);

% Surface temperature [degrees C]
T_s = PM_inp.T_a + PM_inp.H./(rho_a.*c_p.*G_a);

% Calculate Jarvis and McNaughton decoupling coefficient
Omega = (Delta + gma)./(Delta + gma.*(1 + G_a./G_c_pm_v));

% PART 4: INVERT MCCOLL (2020) TO FIND CANOPY CONDUCTANCE, G_c
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

if isequal(PM_inp.invertMcColl,1)
    
    % Convert vapor pressure to specific humidity [kg H20/kg air]
    q_sat =  e_sat.*((M_w/M_a)./PM_inp.P_a);
    q_a = q_sat.*PM_inp.RH;
    
    steps = length(PM_inp.time);
    G_c_mc_v = zeros(steps,1);
    
    for step = 1:steps
        f = @(G_c) invert_McColl(G_c,step);
        if isnan(G_c_pm_v(step)) | G_c_pm_v(step) < 0
            x_i = 0.01;
        else
            x_i = G_c_pm_v(step);
        end
        
        try
            % Numerical solution for canopy/surface conductance from McColl
            % (2020) version of Penman-Monteith [m/s]
            G_c_mc_v(step,1) = fzero(f,x_i);
        catch
            G_c_mc_v(step,1) = NaN;
        end

        display(step/steps*100)
    end
    
    
    % Convert from volume [m/s] to molar units [moles air/m^2 GA/s]
    G_c_mc = MpsFlux2Mol(G_c_mc_v,PM_inp.P_a,PM_inp.T_a);
    
    % Aggregate G_c_pm to daily PET in depth units [mm/day]
    G_c_mc_d = aggData(PM_inp.time,G_c_mc,PM_inp.dt,'mean');
    
    % Calculate an estimate for VPD at the leaf [Pa]
    VPD_l_mc = PM_inp.LE.*gma./(rho_a.*c_p.*G_c_mc_v);

else
    G_c_mc = [];
    G_c_mc_d = [];
    VPD_l_mc = [];
end


    function Res = invert_McColl(G_c,step)
    
        G_w = G_a(step).*G_c./(G_a(step) + G_c);
        
        % Messy Lambert W portion
        W_o = lambertw(((L_v^2*q_sat(step))./(c_p*R_v*(PM_inp.T_a(step) + 273.15).^2).*1./(G_a(step)./G_c + 1)).*...
            exp(L_v./(R_v*(PM_inp.T_a(step) + 273.15).^2).*(Q_ne(step) + G_w.*rho_a(step).*L_v.*q_a(step))./(rho_a(step).*c_p.*G_a(step))));
        
        % McColl (2020) revised Penman-Monteith ET estimate [PET_pt; mm/day]
        LE_mc = (rho_a(step).*c_p.*G_a(step).*W_o./(L_v./(R_v*(PM_inp.T_a(step) + 273.15).^2)) -...
            rho_a(step).*L_v.*G_w.*q_a(step));
        
        % Difference in modelled and observed LE
        Res = LE_mc - PM_inp.LE(step);
        
    end
    

end
