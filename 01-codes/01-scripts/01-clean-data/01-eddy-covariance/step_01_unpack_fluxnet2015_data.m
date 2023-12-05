% =========================================================================
% Name   : step_01_unpack_fluxnet2015_data.m
% Author : Brandon Sloan
% Date   : 7/5/22
%
% DESCRIPTION
% This script upacks the half-hourly/hourly FLUXNET2015 data files to
% nicely formatted timetables in MATLAB.
%   
% =========================================================================
clc
clear
fpath = '.\02-data\01-raw\01-eddy-covariance\01-fluxnet2015\';
folders = dir([fpath,'FLX*']);
spath = '.\02-data\02-processed\01-eddy-covariance\';

for ii = 1:length(folders)
    try
    [Obs_HH,name] = importFluxnetData([fpath,folders(ii).name],'HH');
    catch
    [Obs_HH,name] = importFluxnetData([fpath,folders(ii).name],'HR'); 
    end
    save([spath,name],'Obs_HH')
end
