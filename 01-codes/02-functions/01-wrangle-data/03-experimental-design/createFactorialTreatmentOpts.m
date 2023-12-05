% =========================================================================
% Name   : Create_Tier1_Treatment_Opts.m
% Author : Brandon Sloan
% Date   : 8/2/22
%
% DESCRIPTION
% This script carries creates separate opts files for each "Treatment"
% scenario that I am running for Tier 1.  
% =========================================================================

function [Fit_opts,PP_opts,TrtCode] = createFactorialTreatmentOpts(ED,Fit_opts,PP_opts)
    
% Treatment 1: Gc equation to use
%''''''''''''''''''''''''''''''''
switch ED(1)
    case 1 % Lin
        Fit_opts.GceqnFlag = 1;
    case 2 % Medlyn
        Fit_opts.GceqnFlag = 2;
end

% Treatment 2: Flux2Fit
%''''''''''''''''''''''

switch ED(2)
    case 1 % Gc
        Fit_opts.Flux2Fit = 'G_c';
    case 2 % ET
        Fit_opts.Flux2Fit = 'ET';
end


% Treatment 3: Parameters to fix
%'''''''''''''''''''''''''''''''

switch ED(3)
    case 1 % All
        Fit_opts.prmfixFlag = 1;
    case 2 % G1, m
        Fit_opts.prmfixFlag = 2;
    case 3 % Go, G1
        Fit_opts.prmfixFlag = 3;
    case 4 % G1
        Fit_opts.prmfixFlag = 4;
        
end


% Treatment 4: VPD to fit with
%'''''''''''''''''''''''''''''

switch ED(4)
    case 1 % VPD_a
        Fit_opts.VPDFlag = 0;
    case 2 % VPD_l
        Fit_opts.VPDFlag = 1;
    case 3 % Iterate VPD_l
        Fit_opts.VPDFlag = 2;    
end


% Treatment 5: Fitting settings
%''''''''''''''''''''''''''''''

switch ED(5)
    case 1 % NLS
        Fit_opts.SolverFlag = 1;
        Fit_opts.Robust = 0;
        Fit_opts.WeightFlag = 0;
        Fit_opts.LossFxnFlag = '';
    case 2 % Robust NLS
        Fit_opts.SolverFlag = 1;
        Fit_opts.Robust = 1;
        Fit_opts.WeightFlag = 0;
        Fit_opts.LossFxnFlag = '';
    case 3 % Constrained, L1 Weighted optimization
        Fit_opts.SolverFlag = 4;
        Fit_opts.Robust = 0;
        Fit_opts.WeightFlag = 1;
        Fit_opts.LossFxnFlag = 'L1';  
    case 4 % Constrained, LCE optimization
        Fit_opts.SolverFlag = 4;
        Fit_opts.Robust = 0;
        Fit_opts.WeightFlag = 0;
        Fit_opts.LossFxnFlag = 'LCE';
end


% Treatment 6: SEB correction
%''''''''''''''''''''''''''''''

switch ED(6)
    case 1 % No SEB correction
        PP_opts.SEBFlag = 0;
    case 2 % Bowen ratio
        PP_opts.SEBFlag = 1;
    case 3 % All error in LE
        PP_opts.SEBFlag = 2;
%     case 4 % All error in Rn
%         PP_opts.SEBFlag = 3;
end


% Treatment 7: Growing seasona filter
%''''''''''''''''''''''''''''''''''''

switch ED(7)
    case 1 % No GPP filter
        PP_opts.GPPFiltFlag = 0;
%     case 2 % Remove < 10% of 95th prctile
%         PP_opts.GPPFiltFlag = 1;
    case 2 % Remove < 50% of 95th prctile
        PP_opts.GPPFiltFlag = 2;
end

% Treatment 8: Include LAI
%'''''''''''''''''''''''''
switch ED(8)
    case 1 % No LAI
        PP_opts.LAIFlag = 0;
    case 2 % LAI
        PP_opts.LAIFlag = 1;
end

TrtCode = getCode(Fit_opts,PP_opts);

end

function TrtCode = getCode(Fit_opts,PP_opts)

% Set default name to save
switch Fit_opts.GceqnFlag
    case 1
        Gceqn = 'Lin';
    case 2
        Gceqn = 'Med';
end

TrtCode = [Gceqn,strrep(Fit_opts.Flux2Fit,'_',''),'_P',...
        num2str(Fit_opts.prmfixFlag),'_VPD',num2str(Fit_opts.VPDFlag),...
        'O',num2str(Fit_opts.SolverFlag),'n',num2str(Fit_opts.nboot),...
        '_EB',num2str(PP_opts.SEBFlag),'_GPP',num2str(PP_opts.GPPFiltFlag),...
        '_LAI',num2str(PP_opts.LAIFlag)];

end


