function [Obs,PM_opts] = getBaseFluxDataforFactorial(SiteProp,PP_opts,PM_opts)

% Extract yearly or monthly filtering criteria if required for crops
switch SiteProp.h_vSourceFlag
    case 4 % Year or monthly criteria
        PP_opts.yrs = SiteProp.CropYears{:};
        PP_opts.mns = SiteProp.CropMonths{:};
    case 5 % Select dominant crop of rotation
        if iscell(SiteProp.CropYears{:})
            PP_opts.yrs = SiteProp.CropYears{1,1}{1,1};
        else
            PP_opts.yrs = SiteProp.CropYears{:};
        end
        if iscell(SiteProp.CropMonths{:})
            PP_opts.mns = SiteProp.CropMonths{1,1}{1,1};
        else
            PP_opts.mns = SiteProp.CropMonths{:};
        end
    otherwise
        PP_opts.yrs = [];
        PP_opts.mns = [];
end

% Unpack and filter flux data
PP_opts.sm_sensor = SiteProp.SWC_sensor_used{:};
[Obs,dt,theta_sat,PM_opts.theta_prc,remfilt] = ...
    preprocessFluxData(SiteProp.Site_ID{:},PP_opts);

% Gapfill CO2 data if necessary
if sum(isnan(Obs.C_a))/size(Obs,1) > 0.2
   load NOAA_Mauna_Loa_Annual_CO2.mat
   yrs = year(Obs.Time);
   [co2match,idx] = ismember(yrs,ML_CO2.Year);
   Obs.C_a(co2match) = ML_CO2.CO2_ppm(idx);
else
end

% Allows fit analysis without soil water potential
try
    % Calculate soil water potential
    Obs.psi_s = theta2psi(SiteProp.psi_sat,SiteProp.b,Obs.SWC./theta_sat);
catch
end

% Invert Penman-Moteith to get Gc
PM_opts.dt = dt;
PM_opts.h_v = SiteProp.h_v;
PM_opts.z_m = SiteProp.z_m;
PM_opts.z_0m_c = PM_opts.h_v*0.1;
Obs = calcGcfromFluxdata(Obs,PM_opts);

% Pass GPP days to keep correct filtering from original time series
PM_opts.GPPdays = remfilt.GPPdays;

end