% =========================================================================
% Name   : step_01_fit_pmoc_parameters.m
% Author : Brandon Sloan
% Date   : 6/1/22
%
% DESCRIPTION
% This script fits the parameters of the Penman-Monteith Optimal
% Conductance (PMOC) model under the 2,304 unique data filtering and model
% specification assumptions described in Sloan and Feng (2023). This script
% is a shell for FLUXNET_bootfit_Factorial_Fast.m that does the
% model-fitting and file saving. The code is parallelizable and I have
% provided instructions on how to perform this either in this script or in
% downstream scripts if of interest. The fits can take quite a while per
% site (~10 minutes) when run serially as there are over 2000 nonlinear 
% regressions taking place.
%
% =========================================================================

clc
clear
close all

% Note: Uncomment these lines and insert parfor on line 51 of this file or
% line 32 of fit_pmoc_factorial.m
% % Add the necessary paths if parallel runs are on workers without the path
% % information.
% addpath(genpath('./'))
% 
% Generate the parallel pool. Sometimes the parallel pool will not start,
% so this attempts to run it again if it fails.
workers = 12; parpool(workers)
if  isempty(gcp('nocreate'))
     parpool(workers)
     if isempty(gcp('nocreate'))
         error('Could not get all workers')
     else
     end
else
end

% Load data
load final_ec_site_properties.mat
load pmoc_default_opts.mat

% Set and store initial settings
Fit_opts_o = Fit_opts;
PP_opts_o = PP_opts;
PM_opts_o = PM_opts;

% DO NOT TOUCH. Create experimental design matrix. These options are 
% translated into differing data filtering and model specifications 
% during the fit.
ED = fullfact([2,2,4,3,4,3,2,2]);

% Insert parfor here if desired. Also, this code will only run for the
% given example files US-Me1 and US-Me2
for ii = 126

% Fit PMOC model for all assumption sets
fit_pmoc_factorial(ED,SiteProp(ii,:),PP_opts_o,PM_opts_o,...
    Fit_opts_o,'check');

end
