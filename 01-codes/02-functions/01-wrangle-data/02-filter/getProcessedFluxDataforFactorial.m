function Obs = getProcessedFluxDataforFactorial(Obs,PP_opts,PM_opts)

% Prepocess obs futher for SEB options 3
switch PP_opts.SEBFlag
    
    case 0 % All error in H (as is)
        
    case{1,2,3}
        
        switch PP_opts.SEBFlag
            case 1 % Use given FLUXNET bowen ratio correction
                    Obs.LE = Obs.LE_bc;
                    
            case 2 % All error in LE
                Obs.LE = Obs.R_n - Obs.H - Obs.G;
                
            case 3 % All error is in R_n
                Obs.R_n = Obs.LE + Obs.H + Obs.G;
        end
        
        % Re-calculate Gc and filter data
        Obs = calcGcfromFluxdata(Obs,PM_opts);
end


% Remove any new NaN's that could spring up
nankeep = sum(isnan(Obs{:,PP_opts.nanchkNames}),2) == 0 & Obs.LE > ...
    PP_opts.LE_gt & Obs.R_n > PP_opts.R_n_gt;
Obs = Obs(nankeep,:);

% Prepocess Obs further for GPP
switch PP_opts.GPPFiltFlag
    case 0
    case {1,2}
        hhtime = dateshift(Obs.Time,'start','day');
        GPPkeep = ismember(hhtime,PM_opts.GPPdays);
        Obs = Obs(GPPkeep,:);
end

% Preprocess further for GPP if necessary
switch PP_opts.LAIFlag
    case 0
        LAIkeep = ones(size(Obs,1),1);
    case 1
        Obs{:,PP_opts.LAIvars} = Obs{:,PP_opts.LAIvars}./Obs.LAI;
        LAIkeep = Obs.LAI > 0 & ~isnan(Obs.LAI);
end

% Remove outliers, threshold violations and LAI errors if necessary
Gckeep = ~isoutlier(Obs.G_c,PP_opts.OutlierMethod);
thkeep = ~isnan(Obs.LE) & ...
    Obs.VPD_l > 0 & ~isinf(Obs.VPD_l) & ~isnan(Obs.VPD_l);
Obs = Obs(Gckeep & thkeep & LAIkeep,:);

end