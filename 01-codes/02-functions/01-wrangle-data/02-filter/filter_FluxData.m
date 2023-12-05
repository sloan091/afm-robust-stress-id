% =========================================================================
% Name   : filter_FluxData.m
% Author : Brandon Sloan
% Date   : 7/27/21
%
% DESCRIPTION
% This function removes poor quality and NaN data points from the flux data
% used to fit the Medlyn-like G_c models.  Flux- and state-specific
% thresholds are used to exclude times of minimal plant transpiration and
% CO2 assimilation (e.g., cold periods, nighttime). These criteria are
% based on the suggestions of Li et al. (2018) and Lin et al. (2019).
% =========================================================================

function [Obsfilt,remfilt] = filter_FluxData(PP_opts,Obs)

% Dummy keep logical vector used throughout
dmykeep = true(size(Obs,1),1);

% Remove NaN time steps
remNaNall = isnan(Obs{:,PP_opts.nanchkNames});
[nanlimct,maxid] = max(sum(remNaNall,1));
nanlimvar = PP_opts.nanchkNames{maxid};
nankeep = sum(remNaNall,2) == 0;

% Remove poor quality data with QC flags
QCNames = append(PP_opts.QCFlagNames,'_QC');
remQCall = Obs{:,QCNames} > PP_opts.QCFlagMax & ~isnan(Obs{:,QCNames});
[QClimct,maxid] = max(sum(remQCall,1));
QClimvar = QCNames{maxid};
QCkeep = sum(remQCall,2) == 0;

% Change key variable names
old = Obs.Properties.VariableNames;
[match,midx] = ismember(PP_opts.ObsNames,old);
Obs.Properties.VariableNames(midx(match)) = PP_opts.NewNames(match);

% Select only realistic fluxes and states
threshkeep = Obs.LE > PP_opts.LE_gt  & Obs.GPP > PP_opts.GPP_gt & ...
    Obs.H > PP_opts.H_gt & Obs.R_n > PP_opts.R_n_gt & Obs.RH < PP_opts.RH_lt ...
    & ismember(Obs.NIGHT,0) & Obs.T_a > PP_opts.T_a_gt ...
    & Obs.VPD_a > PP_opts.VPD_gt & sum(isnan(Obs{:,{'LE','GPP','H','R_n',...
    'RH','NIGHT','T_a','VPD_a'}}),2) == 0;

% Remove specified time steps after rainfall
Prem = Obs.(PP_opts.PrecipName) > 0;
Pind = find(Prem);
dur = PP_opts.Prem_dt/PP_opts.dt;
remove = ones(dur,1);
for i = 1:length(Pind)
    Prem(Pind(i):Pind(i)+dur-1,1) = remove;
end
Pkeep = logical(~Prem(1:height(Obs),1));

% Filter out data based on growing season GPP threshold
switch PP_opts.GPPFiltFlag
    case 0 % Special case used to pass the usable days to my factorial code
        PP_optstmp = PP_opts;
        PP_optstmp.GPPFiltFlag = 2;
        [~,~,remfilt.GPPdays] = filtGPP(PP_optstmp,Obs);
         GPPkeep = dmykeep;
    case {1,2}
        [~,GPPkeep,remfilt.GPPdays] = filtGPP(PP_opts,Obs);
end

% Select only specific months and years for crops and/or crop rotations
% Year select
if isempty(PP_opts.yrs)
    yrskeep = dmykeep;
else
    yrskeep = ismember(year(Obs.Time),PP_opts.yrs);
end
% Month select
if isempty(PP_opts.mns)
    mnskeep = dmykeep;
else
    mnskeep = ismember(month(Obs.Time),PP_opts.mns);
end


% Select final data set by removing bad rows and unnecessary columns
keep = nankeep & QCkeep & threshkeep & Pkeep & GPPkeep & yrskeep & mnskeep;
Obsfilt = Obs(keep,midx(midx>0));

% Store removed data counts
remfilt.n_in = size(Obs,1);
remfilt.n_filt1 = size(Obsfilt,1);
remfilt.perc_left = remfilt.n_filt1/remfilt.n_in*100;
remfilt.QClim = QClimct./remfilt.n_in*100;
remfilt.QClimvar = QClimvar;
remfilt.nanlim = nanlimct./remfilt.n_in*100;
remfilt.nanlimvar = nanlimvar;
remfilt.Plim = sum(~Pkeep)./remfilt.n_in*100;
remfilt.threshlim = sum(~threshkeep)./remfilt.n_in*100;
remfilt.GPPlim = sum(~GPPkeep)./remfilt.n_in*100;

end