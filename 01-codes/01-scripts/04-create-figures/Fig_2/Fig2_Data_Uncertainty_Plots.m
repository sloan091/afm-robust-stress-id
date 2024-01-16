% =========================================================================
% Name   : ET_v_SM_Filtered.m
% Author : Brandon Sloan
% Date   : 2/14/23
%
% DESCRIPTION
% This checks each FLUXNET2015/AmeriFlux-FLUXNET dataset for the requisite
% soil moisture and atmospheric forcing measurments needed for assessing
% the plant response to soil and atmospheric dryness at the flux tower
% sites. The script creates a CandidateSites table identify the limiting
% measurement or inadequate sites due to missing measurements.
%
% =========================================================================

clc
clear
close all

% Load data
load Tier1_Analysis_Default_Opts.mat
load Tier1_FLUXNET2015_SiteProp.mat

% Alter default fit parameters
PP_opts.yrs = [];
PP_opts.mns = [];

% Select sites
nsites = size(SiteProp,1);
names = SiteProp.Site_ID;

% Reset bins
Fit_opts.binedges = 0:0.1:1;
Fit_opts.bins = 10;


Sites = cell(nsites,23);
for ii = 126%:nsites
    try
        % Load raw observations
        Obs = rawLoad(names{ii},PP_opts);
        
        % Load filtered observations
        PP_opts.GPPFiltFlag = 1;
        PP_opts.binedges = Fit_opts.binedges;
        [Obs_filt,dt,theta_sat,theta_prc,remfilt] = ...
            preprocessFluxDataCandidate(names{ii},PP_opts);
        remfilt = rmfield(remfilt,'GPPdays');
        Sites(ii,:) = [names{ii},struct2cell(remfilt)'];
        
        % Load filtered observations and remove outliers
        [Obs_filt,PM_opts] = getBaseFluxDataforFactorial(SiteProp(ii,:),...
            PP_opts,PM_opts);
        Gckeep = ~isoutlier(Obs_filt.G_c,PP_opts.OutlierMethod);
        VPDkeep = Obs_filt.VPD_l > 0 & ~isinf(Obs_filt.VPD_l);
        Obs_filt = Obs_filt(Gckeep & VPDkeep,:);
        
    catch ME
        Sites(ii,1) = names(ii);
        Sites(ii,6) = {ME.message};
        Sites(ii,end-4) = {1};
    end
    
end

close all
    ar = 4/3;
    wd = 1*2.54;
    ht = wd/ar;
    fr = 14/12;
    fs = 6;
    ft = 'CMU Sans Serif';
    %ft = 'Helvetica';
   [h,ax] = plotScatter(Obs_filt.SWC,Obs_filt.LE./(Obs_filt.LE + Obs_filt.H),...
       'Soil Moisture Percentile, \theta_p','ET [mm/d]','WW','o',1,'none',...
       Obs_filt.VPD_a,0.5,[],[],fs/fr,ft,[]);
    set(ax,'XTick',[0,0.5,1],'YTick',[0,6,12]);
    h.MarkerFaceAlpha = 0.2;
    ax.LabelFontSizeMultiplier = fr;
    fh = gcf;
    fh.Position(3)/2.54

% Latent heat flux
close all
    ar = 4/3;
    wd = 1*2.54;
    ht = wd/ar;
    fr = 14/12;
    fs = 6;
    ft = 'CMU Sans Serif';
    %ft = 'Helvetica';
   [h,ax] = plotScatter(Obs_filt.SWC_nep,Obs_filt.LE*86400/2.5e6,...
       'Soil Moisture Percentile, \theta_p','ET [mm/d]','WW','o',1,'none',...
       [0.75,0.75,0.75],0.5,[0,1],[0,15],fs/fr,ft,[0.8,0.8,wd,ht]);
    set(ax,'XTick',[0,0.5,1],'YTick',[0,6,12]);
    h.MarkerFaceAlpha = 0.2;
    ax.LabelFontSizeMultiplier = fr;
    fh = gcf;
    fh.Position(3)/2.54
    print('SMET','-dpng','-r1200')
    
    % Canopy conductance, Gc
    close all
    ar = 4/3;
    wd = 1*2.54;
    ht = wd/ar;
    fr = 14/12;
    fs = 6;
    ft = 'CMU Sans Serif';
    %ft = 'Helvetica';
   [h,ax] = plotScatter(Obs_filt.SWC_nep,Obs_filt.G_c,...
       'Soil Moisture Percentile, \theta_p','G_c [nol/m^2/s]','WW','o',1,'none',...
       Obs_filt.VPD_a,0.5,[0,1],[0,0.6],fs/fr,ft,[0.8,0.8,wd,ht]);
    set(ax,'XTick',[0,0.5,1],'YTick',[0,0.3,0.6]);
    h.MarkerFaceAlpha = 0.2;
    ax.LabelFontSizeMultiplier = fr;
    fh = gcf;
    fh.Position(3)/2.54
    cmocean('thermal');
    print('SMGc','-dpng','-r1200')
    
    
    % SEB Closure
    close all
    SEBfit = fitlm(Obs_filt.R_n + Obs_filt.G,Obs_filt.LE + Obs_filt.H);
    [h,ax] = plotScatter(Obs_filt.R_n + Obs_filt.G,Obs_filt.LE + Obs_filt.H,...
        'R_n - G [W/m^2]','H + LE [W/m^2]','WW','o',1,'none',...
        [0.75,0.75,0.75],0.5,[],[],fs/fr,ft,[1,1,wd,ht]);
    h.MarkerFaceAlpha = 0.2;
    ax.LabelFontSizeMultiplier = fr;
    hold on
    plot([0,1e3],[0,1e3],'Color','k','DisplayName','1:1',...
        'LineStyle','--','LineWidth',1)
    hold on
    plot([0,1e3],SEBfit.predict([0;1e3]),'Color','k','DisplayName','Fit',...
        'LineWidth',1)
    gax = gca;
    gax.XTick = [0,1e3];
    gax.YTick = [0,1e3];
    gax.XTickLabel{2} = '10^3';
    gax.YTickLabel{2} = '10^3';
    print('SEB','-dpng','-r1200')
    
    % Random uncertainty
    close all
    ar = 4/3;
    wd = 1*2.54;
    ht = wd/ar;
    h = histogram(Obs_filt.LE_ru*86400/2.5e6,0:0.2:4,...
        'Normalization','probability','FaceColor',[0.75,0.75,0.75],'EdgeColor','none');
    hold on
    plotLine([],[],'ET noise [mm/d]',...
        'Probability','','-',2,...
        0.5,[0,4],[0,0.13],fs/fr,ft,'b',[1,1,wd,ht]);
        gax = gca;
    gax.LabelFontSizeMultiplier = fr;
    gax.XTick = [0,4];
    gax.YTick = [0,0.1];
    print('ETru','-dpng','-r1200')
    
    % Box chart
    close all
    bins = discretize(Obs_filt.LE,prctile(Obs_filt.LE,0:10:100));
    %[bins,X] = discretize(Obs_filt.LE,10);
    bh = boxchart(bins,Obs_filt.LE_ru*86400/2.5e6);
    bh.BoxFaceColor = [0.75,0.75,0.75];
    bh.BoxFaceAlpha = 0;
    bh.WhiskerLineColor = [0.75,0.75,0.75];
    bh.MarkerStyle = 'none';
    bh.BoxWidth = 0.7;
    
    ar = 4/3;
    wd = 1*2.54;
    ht = wd/ar;
    hold on
    plotLine([],[],'ET percentile',...
        'ET noise [mm/d]','','-',2,...
        0.5,[0,10.5],[0,3],fs/fr,ft,'b',[0.6,0.6,wd,ht]);
        gax = gca;
    gax.LabelFontSizeMultiplier = fr;
    gax.XTick = [0,10];
    gax.YTick = [0,3];
    gax.XTickLabel = {'0','1'};
    %print('ETru','-dpng','-r1200')
    print('ETru_box','-dsvg')
    
    
    % Phenology plot
    close all
    ar = 8/3;
    wd = 2*2.54;
    ht = wd/ar;
    load US_Me2_LAI_FLX15.mat
    %      stepseas = calc_Seasonal_Means(LAI_FLX15.Time,...
    %             LAI_FLX15.LAI,15);
    
    % Grab 95th perctile of daily GPP data
    tmp = load([names{ii},'_Daily']);
    D = tmp.Obs_Daily;
    GPPmax = prctile(D.GPP_NT_VUT_REF,95);
    datewdw = [datetime(2017,1,1),datetime(2020,1,1)];
    
    % Smooth with 15 average filter
    GPPsm = smoothdata(D.GPP_NT_VUT_REF,'movmean',days(15),'omitnan',...
        'SamplePoints',D.Time);
    GPP_rem = GPPsm < 0.5*GPPmax;
    yyaxis right
    plotLine(LAI_FLX15.Time,...
        LAI_FLX15.LAI,'Date',...
        'LAI','LAI','-',1,...
        0.5,[],[],fs/fr,ft,[0.8,0.8,0.8],[])
    gax1 = gca;
    gax1.YTick = [0,2];
    hold on
    yyaxis left
    plotLine(D.Time,GPPsm,'Date',...
        'GPP','GPP','-',1,...
        0.5,datewdw,[],fs/fr,ft,'k',[0.8,0.8,wd,ht]);
    hold on
    plot(datewdw,GPPmax*0.5*ones(2,1),'Color','k','LineStyle',':');
    gax2 = gca;
    f = gcf;
    gax = gca;
    %    gax.XTick = [0,1000];
   gax2.YTick = [0,10];
   datetick('x','mm/yy')
   gax.LabelFontSizeMultiplier = fr;
   
%    f.Position(3) = f.Position(3) + 2*gax.TightInset(1);
%    f.Position(4) = f.Position(4) + 2*gax.TightInset(2);
   print('Pheno','-dsvg')


    
  
    
    
    function Obs = rawLoad(name,PP_opts)
    
    % Load raw observations
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
    
    % Calculate SWC exceedance probabilities
    Obs.SWC_nep = 1 - getPrctile(Obs.(sm_sensor));
    
    % Identify which soil moisture sensor to use
    %     sm_sensor = PP_opts.sm_sensor;
    PP_opts.nanchkNames = [PP_opts.nanchkNames,sm_sensor];
    PP_opts.ObsNames = [PP_opts.ObsNames,sm_sensor];
    PP_opts.NewNames = [PP_opts.NewNames,'SWC'];
    
    % Change key variable names
    old = Obs.Properties.VariableNames;
    [match,midx] = ismember(PP_opts.ObsNames,old);
    Obs.Properties.VariableNames(midx(match)) = PP_opts.NewNames(match);
    
    
    end
  
