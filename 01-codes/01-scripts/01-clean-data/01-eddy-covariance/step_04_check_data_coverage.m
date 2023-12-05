% =========================================================================
% Name   : step_04_check_data_coverage.m
% Author : Brandon Sloan
% Date   : 7/5/22
%
% DESCRIPTION
% This checks each FLUXNET2015/AmeriFlux-FLUXNET dataset for the requisite
% soil moisture and atmospheric forcing measurments needed for assessing
% the plant response to soil and atmospheric dryness at the flux tower
% sites. The script creates a SiteCoverage table identify the limiting
% measurement or inadequate sites due to missing measurements.
%
% =========================================================================

clc
clear
close all

% Load data
load pmoc_default_opts
load potential_ec_sites.mat

% Alter default fit parameters
PP_opts.yrs = [];
PP_opts.mns = [];

% Select sites
nsites = size(SiteList,1);
names = SiteList.Site_ID;

% Reset bins
Fit_opts.binedges = 0:0.1:1;
Fit_opts.bins = 10;


Sites = cell(nsites,25);

% NOTE: the code will only work for US-Me1 and US-Me2 since these are the
% included eddy covariance data. I have included the full site coverage
% table for all sites under the file ec_site_coverage.mat
for ii = 163:164
    try
        % Unpack and filter flux data
        PP_opts.GPPFiltFlag = 0;
        PP_opts.binedges = Fit_opts.binedges;
        [Obs,dt,theta_sat,theta_prc,remfilt] = ...
            preprocessFluxDataCandidate(names{ii},PP_opts);
        try
        remfilt = rmfield(remfilt,'GPPdays');
        catch
        end
        Sites(ii,:) = [names{ii},struct2cell(remfilt)'];
        
    catch ME
        Sites(ii,1) = names(ii);
        Sites(ii,6) = {ME.message};
        Sites(ii,end-6) = {1};
    end
    
end

nms = ['SiteID',fieldnames(remfilt)'];
SiteCoverage = cell2table(Sites);
SiteCoverage.Properties.VariableNames = nms;
% Do not want to overwrite the site coverage file included
%save('ec_site_coverage','SiteCoverage');
