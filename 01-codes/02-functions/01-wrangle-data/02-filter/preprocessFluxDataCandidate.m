function [Obs,dt,theta_sat,theta_prc,remfilt] = preprocessFluxDataCandidate(name,PP_opts)

% Unpack flux data
tmp = load(name);
Obs = tmp.('Obs_HH');
vn = Obs.Properties.VariableNames;

% Find any SWC measurements
scol = contains(vn,'SWC') & ~contains(vn,'QC');
smnames = sort(vn(scol));
smprc = sum(isnan(Obs{:,smnames}))/size(Obs,1)*100;
nsensors = length(smnames);

% Select sensor with most data (least NaNs)
if nsensors == 0
    sm_sensor = 'SWC';
elseif nsensors == 1
    sm_sensor = smnames{:};
else
    [~,sel_sm] = min(smprc);
    if smprc(1) - smprc(sel_sm) > 10
        sm_sensor = smnames{sel_sm};
    else
        sm_sensor = smnames{1};
    end
end


% Add surface energy budget correction if necessary

% Assume NaN G values are negligible
Obs.G_F_MDS(isnan(Obs.G_F_MDS)) = 0;

% Calculate surface energy budget imbalance
TF = Obs.LE_F_MDS + Obs.H_F_MDS ;
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
%     sm_sensor = PP_opts.sm_sensor;
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

% Check for QC variables before filtering
chk = ismember(PP_opts.nanchkNames,Obs.Properties.VariableNames);
if sum(chk) < length(chk)
    Obs = [];
    theta_sat = [];
    theta_prc = [];
    %     remfilt.theta_cv = [];
    % Store removed data counts
    remfilt.n_in = [];
    remfilt.n_filt1 = [];
    remfilt.perc_left = [];
    remfilt.QClim = [];
    remfilt.QClimvar = [PP_opts.nanchkNames{~chk},' does not exist'];
    remfilt.nanlim = [];
    remfilt.nanlimvar = [];
    remfilt.Plim = [];
    remfilt.threshlim = [];
    remfilt.GPPlim = [];
    %remfilt.EBslopeAll = SEBfit.Coefficients.Estimate(2);
    remfilt.EBslope = [];
    remfilt.smnames = smnames;
    remfilt.smprc = smprc;
    remfilt.smsensor = sm_sensor;
    remfilt.yr_range = [];
    remfilt.nyrs = [];
    remfilt.GPPaddlim = [];
    remfilt.errFlag = 2;
    remfilt.smcounts = [];
    remfilt.smgoodbins = 0;
    remfilt.yrsbin = [];
    remfilt.yrmnct = [];
    remfilt.expVar = [];
    remfilt.expVarBins = [];
else
    
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
    Obs_o = Obs;
    [Obs,remfilt] = filter_FluxData(PP_opts,Obs);
    
    % Calculate surface energy budget imbalance for selected data
    TF = Obs.LE + Obs.H ;
    AE = Obs.R_n - Obs.G;
    Rnkeep = ~isnan(TF) & ~isnan(AE);
    SEBfitsel = fitlm(AE(Rnkeep),TF(Rnkeep));
    
    % Store
    %remfilt.EBslopeAll = SEBfit.Coefficients.Estimate(2);
    remfilt.EBslope = SEBfitsel.Coefficients.Estimate(2);
    remfilt.smnames = smnames;
    remfilt.smprc = smprc;
    remfilt.smsensor = sm_sensor;
    % Define years of coverage
    years = unique(year(Obs.Time));
    yrmn = min(years);
    yrmx = max(years);
    remfilt.yr_range = [num2str(yrmn),'-',num2str(yrmx)];
    remfilt.nyrs = length(years);
    
    PP_opts.GPPFiltFlag = 2;
    [Obs_o,remfiltGPP] = filter_FluxData(PP_opts,Obs_o);
    remfilt.GPPlim = remfiltGPP.GPPlim;
    remfilt.GPPaddlim = (size(Obs,1) - size(Obs_o,1))/size(Obs,1)*100;
    remfilt.errFlag = 0;
    
    % Some histograms for better decision-making
    [remfilt.smcounts,~,binID] = histcounts(Obs.SWC_nep,PP_opts.binedges);
    remfilt.smgoodbins = sum(remfilt.smcounts > 30);
    remfilt.yrsbin = (yrmn-0.5):(yrmx+0.5);
    remfilt.yrmnct = histcounts2(year(Obs.Time),month(Obs.Time),...
        remfilt.yrsbin,0.5:12.5);
    
    % Calculate the explainable variance
    evfxn = @(noise,signal) 1 - sum(noise.^2,'omitnan')/length(noise)/...
        var(signal,'omitnan');
    G = findgroups(binID);
    remfilt.expVar = evfxn(Obs.LE_ru,Obs.LE);
    remfilt.expVarBins = NaN(size(remfilt.smcounts));
    remfilt.expVarBins(remfilt.smcounts > 0) = ...
        splitapply(@(x) evfxn(x(:,1),x(:,2)),[Obs.LE_ru,Obs.LE], G);
    
    
end 

end