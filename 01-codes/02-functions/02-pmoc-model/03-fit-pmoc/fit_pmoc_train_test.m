% =========================================================================
% Name   : bootfitGcPM2Fluxdata.m
% Author : Brandon Sloan
% Date   : 7/15/22
%
% DESCRIPTION
% This function bootstraps the combined canopy conductanc (G_c)
% Penman-Monteith model to flux tower data.
% =========================================================================

function [prms_trn,gof_trn,gof_tst,info_trn,info_tst] = ...
    fit_pmoc_train_test(In,Fit_opts,PM_opts,PP_opts)

% Calculate invariant Penman-Monteith inputs
PM_opts.LAIFlag = PP_opts.LAIFlag;
In = preparePMvars(In,PM_opts);

% Number of bootstraps (I do not bootstrap in the AFM paper)
nboot = Fit_opts.nboot;

% Pass along stat and parameter dimensions
nstats = 10;

% Determine parameter dimension (no stress parametrization in AFM paper)
switch Fit_opts.StressFlag
    case 0
        nprms = 3;
    case 1
        nprms = 4;
    case {2,3,4,5}
        nprms = 5;
    case 6
        nprms = 6;
end

% Add two to nprms to account for reduction factor and G_v
ntotprms = nprms + 2;

% Perform m of n bootstrap with training and testing
[prms_trn,gof_trn,gof_tst,info_trn,info_tst] = train_fit(nboot,@trainfxn,...
    @testfxn,In,Fit_opts.ParFlag,ntotprms,nstats);


% Training function
    function  [prms_trn,gof_trn,mdl_trn,info_trn] = trainfxn(In)
        try
            % Train Gc-PM model to flux data
            [prms_trn,mdl_trn,respvar_trn,VPD_l] = fit_pmoc(In,Fit_opts);
            
            % Extra step if VPDFlag = 2 (iteration)
            if isequal(Fit_opts.VPDFlag,2)
                In.(Fit_opts.VPD_name) = VPD_l;
            else
            end
            
            % Scale ground area results to per unit leaf area if requested
            if isequal(PP_opts.LAIFlag,1)
                In.(Fit_opts.Flux2Fit) = In.(Fit_opts.Flux2Fit).*In.LAI;
                respvar_trn = respvar_trn.*In.LAI;
                prms_trn(:,1) = prms_trn(:,1).*mean(In.LAI,'omitnan');
                In.GPP  = In.GPP.*In.LAI;
            else
            end
            
            % Calculate reduction factor, G_1/VPD_l^m
            redfac = prms_trn(2)./In.(Fit_opts.VPD_name).^prms_trn(3);
            
            % Calculate vegetation conductance, G_v
            if isequal(Fit_opts.GceqnFlag,1)
                G_v = redfac.*In.GPP./In.C_a;
            else
                G_v = (redfac + 1).*1.6.*In.GPP./In.C_a;
            end
            
            % Store trained parameters/variables
            prms_trn = [prms_trn,median(redfac,'omitnan'),median(G_v,'omitnan')];
            
            % Calculate goodness-of-fit statistics
            [RMSE,NSE,Pbias,R,cRMSE,StdRatio,MeanRatio] = ...
                calc_fit_stats(In.(Fit_opts.Flux2Fit),respvar_trn);
            gof_trn = [RMSE,NSE,Pbias,R,cRMSE,StdRatio,MeanRatio,...
                mean(In.(Fit_opts.Flux2Fit),'omitnan'),mdl_trn.AIC,mdl_trn.BIC];
            
            % Calculate information performance metrics
            info_trn = [1,calcInfoPerformance(In.SWC_nep,In.VPD_a,...
                In.(Fit_opts.Flux2Fit),respvar_trn,10,2),2,...
                calcInfoPerformance(In.SWC_nep,In.GPP,...
                In.(Fit_opts.Flux2Fit),respvar_trn,10,2),3,...
                calcInfoPerformance(In.GPP,In.VPD_a,...
                In.(Fit_opts.Flux2Fit),respvar_trn,10,2)];
        catch
            prms_trn = NaN(1,nprms + 2);
            gof_trn = NaN(1,nstats);
            mdl_trn = [];
            tmp = NaN(1,7);
            info_trn = [1,tmp,2,tmp,3,tmp];
        end
    end

% Testing function
    function  [Out,info_tst] = testfxn(mdl_tst,In)
        try
            
            % Create prediction for test inputs
            In_array = table2array(In(:,mdl_tst.predvars));
            respvar_tst = mdl_tst.fit_mdl(In_array);
            
            % Scale ground area results to per unit leaf area
            if isequal(PP_opts.LAIFlag,1)
                In.(Fit_opts.Flux2Fit) = In.(Fit_opts.Flux2Fit).*In.LAI;
                respvar_tst = respvar_tst.*In.LAI;
            else
            end
            
            % Calculate goodness-of-fit statistics
            [RMSE,NSE,Pbias,R,cRMSE,StdRatio,MeanRatio] = ...
                calc_fit_stats(In.(Fit_opts.Flux2Fit),respvar_tst);
            Out = [RMSE,NSE,Pbias,R,cRMSE,StdRatio,MeanRatio...
                mean(In.(Fit_opts.Flux2Fit),'omitnan'),NaN,NaN];
            
            % Calculate information performance metrics
            info_tst = [1,calcInfoPerformance(In.SWC_nep,In.VPD_a,...
                In.(Fit_opts.Flux2Fit),respvar_tst,10,2),2,...
                calcInfoPerformance(In.SWC_nep,In.GPP,...
                In.(Fit_opts.Flux2Fit),respvar_tst,10,2),3,...
                calcInfoPerformance(In.GPP,In.VPD_a,...
                In.(Fit_opts.Flux2Fit),respvar_tst,10,2)];
        catch
            Out = NaN(1,nstats);
            tmp = NaN(1,7);
            info_tst = [1,tmp,2,tmp,3,tmp];
        end
    end

end

