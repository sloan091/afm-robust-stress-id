% =========================================================================
% Name   : calcMedlynGc.m
% Author : Brandon Sloan
% Date   : 3/25/21
%
% DESCRIPTION
% This function calculates the canopy conductance (G_c) using the
% formulations laid out in either Lin et al. (2018) or Medlyn et al.
% (2017).
%
% INPUTS
%   G_o   - Minimum stomatal conductance for H2O [moles air/m^2/s]
%   G_1   - Empirical Medlyn parameter [kPa^0.5]
%   m     - Exponent of VPD_l [-]
%   VPD_l - Vapor pressure difference between inside the leaf and at
%           the leaf surface [kPa]
%   GPP   - Estimate of net CO2 assimilation [micromoles CO2/m^2/s]
%   C_a   - CO2 concentration [micromoles CO2/moles air]
%
% OUTPUTS
%   G_c   - Canopy conductance of H20 [moles air/m^2/s]
%
% REFERENCES
%   (1) Lin, C., Gentine, P., Huang, Y., Guan, K., Kimm, H., & Zhou, S.
%   (2018). Diel ecosystem conductance response to vapor pressure deficit
%   is suboptimal and independent of soil moisture. Agricultural and
%   Forest Meteorology, 250–251, 24–34.
%   https://doi.org/10.1016/J.AGRFORMET.2017.12.078
%
%   (2) Medlyn, B. E., De Kauwe, M. G., Lin, Y.-S., Knauer, J., Duursma,
%   R. A., Williams, C. A., et al. (2017). How do leaf and ecosystem
%   measures of water-use efficiency compare? New Phytologist, 216(3),
%   758–770. https://doi.org/10.1111/nph.14626
% =========================================================================

function fit_pmoc_factorial(ExpDes,SiteProp,PP_opts_o,PM_opts_o,...
    Fit_opts_o,svname_in)

% Select sites
nexp = size(ExpDes,1);
% Specify number of parameters to fit
ngof = 10;
switch Fit_opts_o.StressFlag
    case 0
        nprms_trn = 5;
    case 1
        nprms_trn = 4;
    case {2,3,4,5}
        nprms_trn = 5;
    case {6,7}
        nprms_trn = 6;
end

% Pre-allocate for storage
prms_trn = NaN(Fit_opts_o.nboot.*Fit_opts_o.bins,nprms_trn,nexp);
gof_trn = NaN(Fit_opts_o.nboot.*Fit_opts_o.bins,ngof,nexp);
gof_tst = gof_trn;
tmpmat = [ones(Fit_opts_o.nboot,1,Fit_opts_o.bins),NaN(Fit_opts_o.nboot,7,Fit_opts_o.bins)];
info_trn = repmat([tmpmat;2*tmpmat;3*tmpmat],1,1,1,nexp);
info_tst = info_trn;
bin_info = cell(nexp,6);

% Load base observations
[Obs_base,PM_opts] = getBaseFluxDataforFactorial(SiteProp,PP_opts_o,PM_opts_o);
PP_opts_o.Site_ID = SiteProp.Site_ID{:};

tic
parfor ii = 1:nexp
    
    [Fit_opts,PP_opts,~] = ...
        createFactorialTreatmentOpts(ExpDes(ii,:),Fit_opts_o,PP_opts_o);
    
    % If weighting the fit, add LE_RANDUNC to the NaN check
    PP_opts.nanchkNames = {'LE','R_n'};
    if ismember(Fit_opts.SolverFlag,[3,4]) && Fit_opts.WeightFlag == 1
        PP_opts.nanchkNames = [PP_opts.nanchkNames,Fit_opts.WeightVar];
    else
    end
    
    % Update settings based on input flag
    switch Fit_opts.VPDFlag
        case 0 % Use VPD_a
            Fit_opts.VPD_name = 'VPD_a';
            Fit_opts.VPDIterFlag = 0;
        case 1 % Use VPD_l
            Fit_opts.VPD_name = 'VPD_l';
            Fit_opts.VPDIterFlag = 0;
        case 2 % Use VPD_l and iterate Gc and PM solutions to balance VPD
            Fit_opts.VPD_name = 'VPD_l';
            Fit_opts.VPDIterFlag = 1;
    end
    
    % Re-process and filter data for specific scenario. For example, GPP
    % filtering will remove more points than not filtering.
    Obs = Obs_base;
    Obs = getProcessedFluxDataforFactorial(Obs,PP_opts,PM_opts);
  
    try
        % Bin fits by variable
        [prms_trn_ii,gof_trn_ii,gof_tst_ii,bin_info_ii,info_trn_ii,info_tst_ii] = ...
            fit_pmoc_sm_bin(Obs,Fit_opts,PM_opts,PP_opts);
        
        % Store results
        prms_trn(:,:,ii) = conv3Dto2D(prms_trn_ii);
        gof_trn(:,:,ii) = conv3Dto2D(gof_trn_ii);
        gof_tst(:,:,ii) = conv3Dto2D(gof_tst_ii);
        info_trn(:,:,:,ii) = info_trn_ii;
        info_tst(:,:,:,ii) = info_tst_ii;
        bin_info(ii,:) = [ii,bin_info_ii];
        
        
    catch ME
        % Throw out warning with site name
        warning([PP_opts_o.Site_ID,' failed due to the following error: ',ME.message])
        
    end
end
toc

% Reformat output files to tables
[prms_trn,gof_trn,gof_tst,info_trn,info_tst] = ...
    reformatFactorialFitData(prms_trn,gof_trn,gof_tst,info_trn,...
    info_tst,1,Fit_opts_o.bins,Fit_opts_o.nboot);

% Add bins to savename
svname = [SiteProp.Site_ID{:},'_',svname_in];
save(svname,'prms_trn','gof_trn','gof_tst','bin_info',...
    'info_trn','info_tst');

end

