% =========================================================================
% Name   : step_02_unpack_ameriflux_fluxnet2015_data.m
% Author : Brandon Sloan
% Date   : 7/5/22
%
% DESCRIPTION
% This script upacks the AmeriFlux FLUXNET product sites. This products
% uses the ONEFLUX processing pipeline for QAQC on AmeriFlux products. This
% code is identical to step_01_unpack_fluxnet2015_data.m, but is included
% because I had originally downloaded the files to different locations.
%   
% =========================================================================
clc
clear
fpath = '.\02-data\01-raw\01-eddy-covariance\02-ameriflux-fluxnet\';
folders = dir([fpath,'AMF*']);
spath = '.\02-data\02-processed\01-eddy-covariance\';

for ii = 1:length(folders)
    try
    [Obs_HH,name] = importFluxnetData([fpath,folders(ii).name],'HH');
    catch
    [Obs_HH,name] = importFluxnetData([fpath,folders(ii).name],'HR'); 
    end
    save([spath,name],'Obs_HH')
end
