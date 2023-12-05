% =========================================================================
% Name   : step_02_filter_modis_lai.m
% Author : Brandon Sloan
% Date   : 6/1/22
%
% DESCRIPTION
% This script takes the raw MODIS LAI data and applies a smoothing
% short-term and long-term filter similar to Ukkola et al., 2021. Here is
% the reference:
%
% Ukkola, A. M., Abramowitz, G., & De Kauwe, M. G. (2022). A flux tower 
% dataset tailored for land model evaluation. Earth Syst. Sci. Data, 14, 
% 449–461. https://doi.org/10.5194/essd-14-449-2022
% =========================================================================

clc
clear

load final_ec_site_properties.mat
svpath = '.\02-data\02-processed\02-modis-lai\';

% FILTERING OPTIONS
% Maximum interpolation gap 
gapmax = duration(24*14,0,0);

% Window for short-term smoothing [days]
swin = 90;

% Window for anomaly smoothing [days]
anwin = [45,45];

% In case there are sites that do not work
badsites = {};
for ii = 125:126
    
    % Load FLUXNET2015
    name = SiteProp.Site_ID{ii};
    
    try
        % Load FLUXNET2015
        name = SiteProp.Site_ID{ii};
        tmp = load(name);
        Obs = tmp.('Obs_HH');
        vn = Obs.Properties.VariableNames;
        
        % Load MODIS LAI
        path = '.\02-data\01-raw\02-modis-lai\';
        LAI = importMODISDataSiteID(path,name,'LAI');
        LAI = LAI(:,5:6);
        
        % Spline interpolate skipping gaps greater than 2 weeks
        LAIsp = retime(LAI,'Daily','fillwithmissing');
        LAIsp = fillmissing(LAIsp,'pchip','MaxGap',gapmax);
        
        % Detect outliers
        snr = min(LAIsp.LAI_mean./LAIsp.LAI_std,10);
        ol = isoutlier(LAIsp.LAI_mean,'movmedian',swin);
        remove = LAIsp.LAI_mean < 0 | ol;
        LAIsp{remove,:} = NaN;
        LAIsp = fillmissing(LAIsp,'pchip','MaxGap',gapmax);
        
        % Smooth data with lowess and 3 month window
       LAIsp.LAIsm = smoothdata(LAIsp.LAI_mean,'rloess',swin);

        
        % Calculate climatology based on day of year
        LAIclim = groupsummary(LAIsp,'Time','dayofyear','mean');
        
        % Replace large gaps with climatology
        gaps = isnan(LAIsp.LAI_mean);
        msdays = day(LAIsp.Time(gaps),'dayofyear');
        doykey = double(LAIclim.dayofyear_Time);
        LAIsp.LAIsm(gaps) = interp1(doykey,LAIclim.mean_LAIsm,msdays);
        
        % Interpolate climatology to time series
        LAIclimts = interp1(doykey,LAIclim.mean_LAIsm,day(LAIsp.Time,'dayofyear'));
        
        % LAI anamoly
        LAIanom = LAIsp.LAIsm - LAIclimts;
        
        % Smooth anomalies with ~3 month moving mean
        LAIanom_sm = smoothdata(LAIanom,'rloess',anwin);
        
        % Add smoothed anomalies back to climatology for finished product
        LAIfinal = LAIclimts + LAIanom_sm;
        
        % Add LAI to FLUXNET2015
        Obs.LAI = interp1(LAIsp.Time,LAIsp.LAIsm,Obs.Time,'nearest');
        
        if Obs.Time(1) < LAIsp.Time(1)
            doykey1 = day(Obs.Time,'dayofyear');
            fix = Obs.Time < LAIsp.Time(1);
            Obs.LAI(fix) = interp1(doykey,LAIclim.mean_LAIsm,doykey1(fix));
        else
        end

        % Internally fill LAI
        plim = @(x) [min(x), max(x)];
        clf
        scatter(LAI.Time,LAI.LAI_mean,'DisplayName','Raw MODIS','Marker','o','SizeData',10,'MarkerEdgeColor','k','MarkerEdgeColor','k')
        hold on
        pkeep = ~ismember(LAIsp.Time,LAI.Time);
        scatter(LAIsp.Time(~ol & pkeep),LAIsp.LAI_mean(~ol & pkeep),'DisplayName','Int. MODIS','Marker','x','SizeData',10,'MarkerEdgeColor','k')
        hold on
        scatter(LAIsp.Time(ol & pkeep),LAIsp.LAI_mean(ol & pkeep),'DisplayName','Outliers','Marker','x','SizeData',10,'MarkerEdgeColor','r')
        hold on
        plot(LAIsp.Time,LAIclimts,'DisplayName','Climatology','LineWidth',2,'Color','b','LineStyle',':')
        hold on
        h = plotLine(Obs.Time,Obs.LAI,'Date',...
            'LAI','Filtered','-',2,...
            1,plim(Obs.Time),[],12,'Arial','r',[]);
       
       % Save figure
       print([svpath,name,'_LAIfilt'],'-djpeg')
       
       % Save ouputs
       LAI_FLX15 = Obs(:,'LAI');
       save([svpath,name,'_LAI_FLX15'],'LAI_FLX15');
       
    catch
        badsites = [badsites;name];
    end
end

