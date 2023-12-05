% =========================================================================
% Name   : step_03_seb_correction.m
% Author : Brandon Sloan
% Date   : 4/27/23
%
% DESCRIPTION
% This script calculates the surface energy budget correction using the
% methods in the OnePipe processing pipeline (EB_CF method 1).  This is
% necessary because some flux sites do not include bias-corrected estimates
% (e.g., LE_CORR), likely because the G flux measurements are missing.
%
% =========================================================================
clc
clear
path = '.\02-data\02-processed\01-eddy-covariance\';
files = dir([path,'*.mat']);
missSEB = {};
tic
for ii = 1:length(files)
    try
        % Unpack flux data
        tmp = load(files(ii).name);
        Obs_HH = tmp.('Obs_HH');
        keep = Obs_HH.LE_F_MDS_QC < 2;
        X = Obs_HH(keep,:);
        
        if isequal(sum(~isnan(X.LE_CORR)),0)
            
            chk = sum(~isnan(X.G_F_MDS));
            
            % Assume NaN G values are negligible
            X.G_F_MDS(isnan(X.G_F_MDS)) = 0;
            
            % Calculate initial correction factors and remove outliers
            TF = X.LE_F_MDS + X.H_F_MDS;
            AE = X.NETRAD - X.G_F_MDS;
            CF = AE./TF;
            rem = isoutlier(CF,'quartiles');
            CF(rem & CF < 0) = NaN;
            
            % Identify hourly period for method 1
            bad_hrs = ismember(hour(X.Time),[3:9,15:21]);
            CF(bad_hrs) = NaN;
            
            % 15 day moving median correction factor approach
            CF_final = movmedian(CF,days(15),'omitnan','SamplePoints',X.Time);
            LE_CORR = X.LE_F_MDS.*CF_final;
            H_CORR = X.H_F_MDS.*CF_final;
            
            % Add corrected data back to original data and save
            Obs_HH.LE_CORR(keep) = LE_CORR;
            Obs_HH.H_CORR(keep) = H_CORR;
            save([path,files(ii).name],'Obs_HH');
            missSEB = [missSEB;{files(ii).name,chk}];
            
        else
        end
        
    catch
        ii
    end
end
toc
% Included with all files
% save('Missing_SEB_Correction','missSEB');