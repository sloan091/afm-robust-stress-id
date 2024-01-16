% =========================================================================
% Name   : Tier1_FLUXNET_Analysis_Create_Default_Opts.m
% Author : Brandon Sloan
% Date   : 7/13/22
%
% DESCRIPTION
% Creates the default pre-processing, Penman-Monteith, and fitting options
% used in my Tier 1 analysis of inferring water stress response from
% FLUXNET sites using optimal canopy conductance (Gc) models coupled with
% Penman-Monteith.
% =========================================================================

clc
clear

% 1. DATA PRE-PROCESSING (PP)
%''''''''''''''''''''''''''''

% Fluxes used as NaN filters
PP_opts.nanchkNames = {'LE_F_MDS','NETRAD','PA_F','TA_F','RH','WS_F',...
    'USTAR','GPP_NT_VUT_REF'};

% Fluxes used for QC flag filters
PP_opts.QCFlagNames  = {'LE_F_MDS','NEE_VUT_REF','TA_F','WS_F'};

% Original flux names to grab from overall flux dataset
PP_opts.ObsNames = {'LE_F_MDS','GPP_NT_VUT_REF','RECO_NT_VUT_REF',...
    'NEE_VUT_REF_RANDUNC','NETRAD','G_F_MDS',...
    'H_F_MDS','PA_F','TA_F','RH','WS_F','USTAR','VPD_F_MDS',...
    'LE_RANDUNC','LE_CORR','LAI','CO2_F_MDS',...
    'SWC_nep','CWD'};

% New names for selected functions used in my analysis
PP_opts.NewNames = {'LE','GPP','R_eco','NEE_ru','R_n','G',...
    'H','P_a','T_a','RH','U','u_star','VPD_a',...
    'LE_ru','LE_bc','LAI','C_a','SWC_nep','CWD'};

% Keep > (gt) or < (lt) of above named fluxes
PP_opts.LE_gt = 0;
PP_opts.GPP_gt = 0;
PP_opts.H_gt = 0;
PP_opts.R_n_gt = 0;
PP_opts.RH_lt = 95;
PP_opts.VPD_gt = 0.5;
PP_opts.T_a_gt = 5;

% Maximum QC flag value: 0. measured only, 1. High quality gap-fill, 2.
% Medium quality gapfill or filled by ERA interim, 3. Poor quality gap-fill
PP_opts.QCFlagMax = 1;

% Hours after rainfall to remove from analysis
PP_opts.Prem_dt = 48;

% Precipitation variable name
PP_opts.PrecipName = 'P_F';

% MATLAB isoutlier detection method:
% 'median','mean','quartiles','grubbs','gesd'
PP_opts.OutlierMethod = 'gesd';

% Placeholder variable used by pre-processing scripts
PP_opts.dt = NaN;

% Non-exceedance probability of soil water content defining saturation
% [0,100]
PP_opts.theta_sat_nep = 95;

% Convert fluxes from per unit ground area to per unit leaf area using
% MODIS leaf area index (LAI) estimates: 0. No, 1. Yes 
PP_opts.LAIFlag = 0;

% Variables to convert using LAI if LAIFlag = 1
PP_opts.LAIvars = {'LE','LE_ru','LE_bc','ET','ET_ru','G_c','GPP',...
    'R_eco','NEE_ru','R_n','G','H'};

% Surface energy budget corrections: 0. No correction (all error in H), 
% 1. Use bowen ratio correction provided by FLUXNET, 2. Add residual to LE,
% 3. Trim R_n to equal H + LE + G.
PP_opts.SEBFlag = 0;

% GPP growing season filter flag: 0. No filter, 1. Remove days less than
% 10% of the 95th percentile of daily GPP as in Lin et al. (2018), 2. Use
% 50% threshold as in Nie et al. (2021)
PP_opts.GPPFiltFlag = 0;


% 2. PENMAN-MONTEITH (PM)
%''''''''''''''''''''''''

% Atmospheric conductance equation (G_a), 1. Traditional log-law, 2.
% Knauer et al. (2018) empirical relationship with stability using u*
PM_opts.CondMethod = 2;

% If CondMethod = 1, specify whether to account for atmopsheric stability
% effects: 0. Neutral conditions, 1. Use Businger-Dyer stability relations
PM_opts.calc_Stability = 0;

% Next 3 options are not used but needed by the invert_Penman_Monteith.m
% function used to grab G_c values.  Keep these set as is
PM_opts.alpha = 1.26;
PM_opts.AggMethod = 'mean';
PM_opts.invertMcColl = 0;

% Placeholder variable used by pre-processing scripts
PM_opts.dt = NaN;
PM_opts.h_v = NaN;
PM_opts.z_m = NaN;
PM_opts.z_0m_c = NaN;


% 3. FIT
%'''''''

% Canopy conductance (G_c) equation options
%''''''''''''''''''''''''''''''''''''''''''

% Flux to fit: 'ET' or 'G_c'
Fit_opts.Flux2Fit = 'G_c';

% Vapor pressure deficit option: 0. VPD_a, 1. VPD_l, 2. Balance VPD_l
% between Gc and PM equations using fixed point iteration.
Fit_opts.VPDFlag = 1;

% If VPDFlag = 2, how many fixed point iterations
Fit_opts.VPDIterN = 5;

% Soil moisture variable: 'SWC_nep', 'SWC', or 'psi_s'
Fit_opts.SM_name = 'SWC';

% G_c equation: 1. Lin (2018), 2. Medlyn (2011,2017)
Fit_opts.GceqnFlag = 1;

% Water stress functions applied to G_1 in G_c equation
% 0. Well-watered, 1. Static exp, 2. Dynamic exp separate VPD, 3. Dynamic
% exp VPD interaction, 4. Static weibull, 5. Dynamic Weibull
Fit_opts.StressFlag = 0;


% Gc parameter settings
%''''''''''''''''''''''

% Which Gc eqn parameters to fit
% 1. Fit all, 2. fix Go, 3. fix m, 4. fix Go,m
Fit_opts.prmfixFlag = 1;

% Soil moisture dependence of G_o if prmfixFlag = 1|3
% 0. Constant value, 1. Linearly declines with SWC_nep
Fit_opts.Go_SMdecline = 0;

% Set G_o for prmfixFlag = 2|4
Fit_opts.G_o = 0;

% Set m for prmfixFlag = 3|4
Fit_opts.m = 0.5;

% Gc parameter initial value bounds
% Order: G_o, G_1, m, a(1), a(2)
Fit_opts.x_li = [0,1,0.5,0,0];
Fit_opts.x_ui = [0.1,4,1,0.1,0.1];

% Gc parameter bounds for constrained fit (same order as initial values)
Fit_opts.x_lb = [0,0,0.01,0,0];
Fit_opts.x_ub = [1,20,1.5,1,1];


% Binning options
%''''''''''''''''
% Variable to bin by
Fit_opts.binvar = 'SWC_nep';

% Bin type: 'percentile' or 'absolute'
Fit_opts.bintype = 'absolute';

% If bin_type = 'absolute', provide the bin edges
Fit_opts.binedges = (0:100/10:100)/100;

% Number of bins if edges or percentiles are not defined
Fit_opts.bins = 10;

% Bin percentiles if bin type is 'percentile'
Fit_opts.binperc = [];


% Bootstrapping
%''''''''''''''

% Parallelize bootstrap: 0. No, 1. Yes
Fit_opts.ParFlag = 0;

% Number of bootstrap replicates
Fit_opts.nboot = 1;


% Solver options
%'''''''''''''''

% Select MATLAB solver
% 1. fitnlm, 2. lsqcurvefit, 3. fminsearch, 4. fmincon
Fit_opts.SolverFlag = 1;

% Use robust fit for SolverFlag = 1 (fitnlm only)
Fit_opts.Robust = 0;

% Default MATLAB robust weighting functions (e.g., bisquare,andrews,etc)
Fit_opts.RobustWgtFun = 'bisquare';

% Perform weighted optimization with SolverFlag = 1,3, and 4.
% Options: 0. Unweighted, 1. 1/WeightVar,  2. 1/WeightVar^2
Fit_opts.WeightFlag = 0;

% Variable to weight fit by
Fit_opts.WeightVar = 'ET_ru';

% If SolverFlag = 3|4, selects the lossfxn: 'L1','L2','LCE'
Fit_opts.LossFxnFlag = 'L1';

% Save options
%'''''''''''''
save('Tier1_Analysis_Default_Opts','PP_opts','PM_opts','Fit_opts');