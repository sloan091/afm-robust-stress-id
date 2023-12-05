function [Obs,dt,theta_sat,theta_prc,remfilt] = preprocessFluxData(name,PP_opts)

    % Unpack flux data
    tmp = load(name);
    Obs = tmp.('Obs_HH');
    
    % Add surface energy budget correction if necessary
    
    % Assume NaN G values are negligible
    Obs.G_F_MDS(isnan(Obs.G_F_MDS)) = 0;
    
    % Calculate surface energy budget imbalance
    TF = Obs.LE_F_MDS + Obs.H_F_MDS;
    AE = Obs.NETRAD - Obs.G_F_MDS;
    Rnkeep = ~isnan(TF) & ~isnan(AE);
    SEBfit = fitlm(AE(Rnkeep),TF(Rnkeep));
    
    switch PP_opts.SEBFlag
        
        case 0 % All error in H (as is)
            
        case 1 % Use given FLUXNET bowen ratio correction
                Obs.LE_F_MDS = Obs.LE_CORR;
            
        case 2 % All error in LE
            Obs.LE_F_MDS = Obs.NETRAD - Obs.H_F_MDS - Obs.G_F_MDS;
            
        case 3 % All error is in R_n
            Obs.NETRAD = Obs.LE_F_MDS + Obs.H_F_MDS + Obs.G_F_MDS;
            
    end
    
    % Identify which soil moisture sensor to use
    sm_sensor = PP_opts.sm_sensor;
    PP_opts.nanchkNames = [PP_opts.nanchkNames,sm_sensor];
    PP_opts.ObsNames = [PP_opts.ObsNames,sm_sensor];
    PP_opts.NewNames = [PP_opts.NewNames,'SWC'];
   
    
    % Load LAI if available
    try
        % Unpack and interpolate phenology data
        tmpname = [name,'_LAI_FLX15'];
        tmp = load(tmpname);
        tmpname = fieldnames(tmp);
        LAI = tmp.(tmpname{:});
        Obs.LAI = LAI.LAI;
    catch
    end
    
    % Site-specific pre-processing inputs
    dt = round((datenum(Obs.Time(2)) - datenum(Obs.Time(1)))*24,...
        1,'significant');
    PP_opts.dt = dt;
    
    % Calculate cumulative water deficit [mm]
    Obs.CWD = calcCWD(Obs.LE_F_MDS*PP_opts.dt*3600/2.5e6,Obs.P_F);
    
    % Calculate SWC exceedance probabilities
    Obs.SWC_nep = 1 - getPrctile(Obs.(sm_sensor));

    % Grab theta_sat
    theta_sat = prctile(Obs.(sm_sensor),PP_opts.theta_sat_nep);
    
    % Grab theta percentiles before filtering
    theta_prc = prctile(Obs.(sm_sensor),0:100);
    
    % Add sitename to PP_opts
    PP_opts.Site_ID = name;
    
    % Pre-process flux data
    [Obs,remfilt] = filter_FluxData(PP_opts,Obs);
    remfilt.EBslope = SEBfit.Coefficients.Estimate(2);
    

end