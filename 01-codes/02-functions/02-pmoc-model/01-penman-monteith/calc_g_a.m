% =========================================================================
% Script Name   : Atmos_Cond.m
% Author        : Brandon Sloan
% Start Date    : Mar 19, 2019
% Last Updated  : Mar 19, 2019
%
% Description   : This function calculates the conductances from the canopy air
% to the measurement height using Monin Obukhov Simlarity Theory (MOST), following the
% common equations used shown by Brutsaert (1982)
%
%   INPUTS:
%   h_v    - Vegetation height [m]
%   L_g    - Obukhov length initial guess; - value indicates unstable, + indicates stable [m]
%   z      - Measurement height [m]
%   z_0m_c - Momentum roughness length for corn canopy from Moneith (2013) [m];
%   LAI    - Leaf area index [m^2 leaf/m^2 ground]
%   SAI    - Stem area index [m^2 stem/m^2 ground] 
%   U      - Streamwise velocity at measurement height [m/s]
%
%   OUTPUTS:
%   g_am - Momentum conductance from the canopy air space to measurement height [m/s]
%   g_ah - Heat conductance from the canopy air space to measurement height [m/s]
%   g_av - Vapor conductance from the canopy air space to measurement height [m/s]
% =========================================================================

function [g_am,g_ah,g_av] = calc_g_a(h_v,z,z_om,U,calc_Stability,T_a,H,u_star)               

% Gravitational constant [m/s^2]
g = 9.8;

% von Karmen constant
k = 0.4; 

% Average density of air (kg/m^3)
rho_a = 1.2;   

% Specific heat at constant pressure [J/kg/K]
c_p = 1004;    

% Heat and vapor roughness length [m]
z_oh = 0.1*z_om;   

% Zero-plane displacement height [m]
d_o = 2/3*h_v;



if isequal(calc_Stability,1)
    
    % Obukhov lenght [m]
    L_obkv = (-T_a.*u_star.^3)./(k*g*(H/(c_p*rho_a)));
    
    % Dimensionless stability parameter
    z_L = (z-d_o)./L_obkv;                               
    
% Stability correction factors taken from Brutsaert (1982)
psi_m = zeros(length(z_L),1);
psi_h = psi_m;
ust = z_L < -0.1;
st = z_L > 0.1;
x1 = (1-16*z_L(ust)).^(0.25);
psi_m(ust) = 2*log((1+x1)/2) + log((1+x1.^2)/2) - 2*atan(x1) + pi()/2;
psi_h(ust) = 2*log((1+x1.^2)/2);
psi_m(st) = -5*z_L(st);
psi_h(st) = -5*z_L(st);

% Resistance formulations from Thom (1972) [s/m]
r_am = (log((z-d_o)/z_om) - psi_m).^2./(U*k^2);
r_ah = (log((z-d_o)/z_om) - psi_m).*(log((z-d_o)/z_oh) - psi_h)./(U*k^2);

else
    
 % Resistance formulations from Thom (1972) [s/m]
r_am = (log((z-d_o)/z_om)).^2./(U*k^2);
r_ah = (log((z-d_o)/z_om)).*(log((z-d_o)/z_oh))./(U*k^2);
   
end

% Convert to conductance form [m/s]
g_am = 1./r_am;
g_ah = 1./r_ah;
g_av = g_ah;

end