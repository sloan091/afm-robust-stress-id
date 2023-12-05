function Obs = calcGcfromFluxdata(Obs,PM_opts)

    % Invert Penman-Monteith
    PM_opts = prepare_PM_inputs(PM_opts,Obs);
    [~,~,G_c_pm,~,VPD_l_pm,~,~,~,T_s,Omega] = ...
        invert_Penman_Monteith(PM_opts);
    
    % Store P-M solution
    Obs.G_c = G_c_pm;
    Obs.VPD_l = VPD_l_pm/1e3;
    
    % Calculate ET from LE [mm/d]
    Obs.ET = Obs.LE*86400/2.5e6;
    
    % Calculate ET random uncertainty [mm/d]
    Obs.ET_ru = Obs.LE_ru*86400/2.5e6;
     
    % Change VPD units to kPa
    Obs.VPD_a = Obs.VPD_a/10;
    
    % Surface temperature
    Obs.T_s = T_s;
    
    % Jarvis and McNaughton decoupling coefficient
    Obs.Omega = Omega;

end