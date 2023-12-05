% =========================================================================
% Name   : binbootfitGcPM2Fluxdata.m
% Author : Brandon Sloan
% Date   : 7/15/22
%
% DESCRIPTION
% This script separates flux tower data into bins and fits the combined
% canopy conductance (G_c) and Penman-Monteith (PM) model to G_c or ET
% observations.
% =========================================================================

function [prms_trn,gof_trn,gof_tst,bin_info,info_trn,info_tst] = ...
    fit_pmoc_sm_bin(In,Fit_opts,PM_opts,PP_opts)

% Create fitting bins
binVar = In.(Fit_opts.binvar);
switch Fit_opts.bintype
    case 'percentile'
        if isempty(Fit_opts.binperc)
            edges = prctile(binVar,[0:100/(Fit_opts.bins):100]);
        else
            edges = prctile(binVar,Fit_opts.binperc);
            Fit_opts.bins = length(Fit_opts.binperc) - 1;
        end
    case 'absolute'
        if isempty(Fit_opts.binedges)
            edges = min(binVar):(max(binVar) - min(binVar))/...
                (Fit_opts.bins):max(binVar);
        else
            edges = Fit_opts.binedges;
            Fit_opts.bins = length(edges) - 1;
        end
end

% LEFT OFF HERE, NEED TO FIX THIS SO I CAN BIN BY YEAR
Fit_opts.TimeBinFlag = 0;
[N,~,idx] = histcounts(binVar,edges);

% A little cheat to get my SWC values when percentiles are used as bins
if strmatch(Fit_opts.binvar,'SWC_nep')
    SWCbins = interp1(0:100,PM_opts.theta_prc,edges*100);
else
    SWCbins = [];
end
bin_info = {Fit_opts.binvar,Fit_opts.bintype,edges,N,SWCbins};

% Pre-allocate
ngof = 10;
switch Fit_opts.StressFlag
    case 0
        nprms_trn = 5;
    case 1
        nprms_trn = 4;
    case {2,3,4,5}
        nprms_trn = 5;
    case {6,7}
        nprms_trn = 6;
end
prms_trn = NaN(Fit_opts.nboot,nprms_trn,Fit_opts.bins);
gof_trn = NaN(Fit_opts.nboot,ngof,Fit_opts.bins);
gof_tst = gof_trn;
tmpmat = [ones(Fit_opts.nboot,1,Fit_opts.bins),NaN(Fit_opts.nboot,7,Fit_opts.bins)];
info_trn = [tmpmat;2*tmpmat;3*tmpmat];
info_tst = info_trn;

% Fit to each bin and store
for ii = 1:Fit_opts.bins
    if N(ii) > 20
        % Select binned data
        select = idx == ii;
        In_sel = In(select,:);
        
        try
            % Bootstrap
            [prms_trn(:,:,ii),gof_trn(:,:,ii),gof_tst(:,:,ii),info_trn(:,:,ii),info_tst(:,:,ii)] =...
                fit_pmoc_train_test(In_sel,Fit_opts,PM_opts,PP_opts);
        catch ME
        end
    else
    end
end

end