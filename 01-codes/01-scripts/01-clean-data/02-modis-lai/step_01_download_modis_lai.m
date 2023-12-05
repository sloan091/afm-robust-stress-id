% =========================================================================
% Name   : step_01_download_modis_lai.m
% Author : Brandon Sloan
% Date   : 5/31/22
%
% DESCRIPTION
% This script automatically downloads the MODIS LAI and fPAR data product
% (MCD15A3H) from the ORNL repository at each eddy covariance site. As a
% bonus I also download NDVI and EVI data products (MOD13Q1), which are not
% used in the AFM paper. Note, to get all sites you must insert the network
% AmeriFlux separately after running for FLUXNET as some of the
% AmeriFlux-FLUXNET sites will not be downloaded.
%   
% =========================================================================
clc
clear

load final_ec_site_properties.mat
q = '"';
base = ['curl -X GET --header ',q,'Accept:text/csv',q,' ',q,'https://modis.ornl.gov/rst/api/v1/'];
prod1 = 'MCD15A3H';
prod2 = 'MOD13Q1';
network = 'FLUXNET/'; % Insert AmeriFlux for the AmeriFlux-FLUXNET sites
dataset = 'subsetStatistics';
daterange = ['?startDate=A1900000&endDate=A2100000',q];
out = ' --output ';
savepath = '.\02-data\01-raw\02-modis-lai\';
tic
% Note only running for the example sites US-Me1 and US-Me2
for ii = 125:126
      
% Load FLUXNET2015
name = SiteProp.Site_ID{ii}; 

curlcmd1 = [base,prod1,'/',network,strrep(name,'_','-'),'/',dataset,daterange,...
    out,savepath,name,'_statistics_',prod1,'.csv'];
curlcmd2 = [base,prod2,'/',network,strrep(name,'_','-'),'/',dataset,daterange,...
    out,savepath,name,'_statistics_',prod2,'.csv'];
[A,~] = system(curlcmd1);
[A,~] = system(curlcmd2);
 
end

dir("*.csv")
toc
