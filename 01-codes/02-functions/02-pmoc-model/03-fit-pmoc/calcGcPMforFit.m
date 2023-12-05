function [estimate,ET,G_c,VPD_l] = calcGcPMforFit(G_o,G_1,m,a,In,Fit_opts)

% Transform back to table
In = array2table(In,'VariableNames',Fit_opts.predvars);

% If Go declines with soil moisture potential
if ismember(Fit_opts.prmfixFlag,[1,3]) && isequal(Fit_opts.Go_SMdecline,1)
    G_o = G_o*(1 - In.SWC_nep);
else
end

% Define forcings used by all settings
VPD = In.(Fit_opts.VPD_name);
GPP = In.GPP;
C_a = In.C_a;
Delta = In.Delta;
gma = In.gma;
Q_ne = In.Q_ne;
ET_i = In.ET_i;
G_a = In.G_a;
c1 = In.c1;

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
rho_a = In.P_a./(R_g*(In.T_a + 273.15))*M_a;


% Soil water stress forcings
if Fit_opts.StressFlag > 0
    SM = In.(Fit_opts.SM_name);
else
    SM = [];
    a = [];
end

% Calculate G_c and set negative values to 0
G_c = calcGc_OptimalWUE(G_o,G_1,m,a,VPD,...
   GPP,C_a,SM,Fit_opts.GceqnFlag,Fit_opts.StressFlag);
G_c(G_c <= 0) = 0;

% Calculate LE with Penman-Monteith convert from W/m^2 to mm/day
LE = (Delta.*Q_ne + ET_i)./...
    (Delta + gma.*(1 + G_a./(G_c.*c1)));
ET = LE*conv;
ET_o = ET;

% Back-calculate VPD_l
VPD_l = LE.*gma./(rho_a.*c_p.*G_c.*c1)/1e3;

% Iterate Gc and ET equations until VPDl converges
if isequal(Fit_opts.VPDIterFlag,1)
        
    for ii = 1:Fit_opts.VPDIterN
        G_c = calcGc_OptimalWUE(G_o,G_1,m,a,VPD_l,...
            GPP,C_a,SM,Fit_opts.GceqnFlag,Fit_opts.StressFlag);
        remGc = G_c <= 0;
        G_c(remGc) = NaN;
        VPD_l(remGc) = NaN;
        LE(remGc) = NaN;
        LE = (Delta.*Q_ne + ET_i)./...
            (Delta + gma.*(1 + G_a./(G_c.*c1)));
        ET = LE*conv;
        % Back-calculate VPD_l
        VPD_l = LE.*gma./(rho_a.*c_p.*G_c.*c1)/1e3;
        
    end
    
    % Set all NaN values to 0
    zeroGc = isnan(G_c);
    G_c(zeroGc) = 0;
    VPD_l(zeroGc) = 0;
    ET(zeroGc) = 0;
    
else
end

switch Fit_opts.Flux2Fit
    case 'ET'
        estimate = ET;
    case 'G_c'
        estimate = G_c;
end
        

end