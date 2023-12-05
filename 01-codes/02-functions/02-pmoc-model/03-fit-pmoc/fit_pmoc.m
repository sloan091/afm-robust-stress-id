% =========================================================================
% Name   : fitGcPM2Fluxdata.m
% Author : Brandon Sloan
% Date   : 7/13/22
%
% DESCRIPTION
% This function fits either Gc or ET to flux tower derived values using the
% combination of ecosystem scale canopy conductance models (G_c) and
% Penman-Monteith for evapotranspriation (ET).  This function is used
% heavily in assessing signals of water stress from FLUXNET data.
% =========================================================================

function [fit_prms,fit_mdl,fit_est,VPD_l] = fit_pmoc(In,Fit_opts)

% 1. SETUP MODEL, VARIABLES, AND PARAMETERS
%''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
% Random start initial x between bounds
rndxi = @(lw,up) (up-lw).*rand(1,length(up)) + lw;
x_i = rndxi(Fit_opts.x_li,Fit_opts.x_ui);

% Define predictor and response variables
predvars = {Fit_opts.VPD_name,'GPP',...
    'C_a','Delta','gma','Q_ne','ET_i','G_a','c1','P_a','T_a'};
respvar = Fit_opts.Flux2Fit;

if Fit_opts.StressFlag > 0
    
    % Add soil moisture to predictors if stressed calculation
    predvars = [predvars,Fit_opts.SM_name];
    
    % Add data if parametrizing soil evaporation under stress
    if isequal(Fit_opts.Go_SMdecline,1)
        predvars = [predvars,'SWC_nep'];
    else
    end
else
end

% Important add predictor names to Fit_opts for later use
Fit_opts.predvars = predvars;

% Set up anonymous estimate functions
pred_fxn = @(G_o,G_1,m,a,In) calcGcPMforFit(G_o,G_1,m,a,In,Fit_opts);


% Determine parameter dimension
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
prmidx = true(1,nprms);
fit_prms = NaN(1,nprms);

% Setup models for solvers
switch Fit_opts.prmfixFlag
    case 1 % Fit Go, G1, m, a
        mdl = @(b,x) pred_fxn(b(1),b(2),b(3),b(4:end),x);
        
    case 2 % Fit G1 and m
        G_o = Fit_opts.G_o;
        mdl = @(b,x) pred_fxn(G_o,b(1),b(2),b(3:end),x);
        prmidx(1) = false;
        fit_prms(1) = mean(G_o,'omitnan');
        
    case 3 % Fit Go and G1
        m = Fit_opts.m;
        mdl = @(b,x) pred_fxn(b(1),b(2),m,b(3:end),x);
        prmidx(3) = false;
        fit_prms(3) = m;
        
    case 4 % Fit G1 only
        m = Fit_opts.m;
        G_o = Fit_opts.G_o;
        mdl = @(b,x) pred_fxn(G_o,b(1),m,b(2:end),x);
        prmidx([1,3]) = false;
        fit_prms([1,3]) = [G_o,m];
        
    case 5 % Fix Go,G1,m
        m = Fit_opts.m;
        G_o = Fit_opts.G_o;
        G_1 = Fit_opts.G_1;
        mdl = @(b,x) pred_fxn(G_o,G_1,m,b(1:end),x);
        prmidx(1:3) = false;
        fit_prms(1:3) = [mean(G_o,'omitnan'),G_1,m];
        
    case 6 % Fix Go,m,exponent
        m = Fit_opts.m;
        G_o = Fit_opts.G_o;
        a = Fit_opts.a;
        mdl = @(b,x) pred_fxn(G_o,b(1),m,[b(2),a],x);
        prmidx([1,3,5]) = false;
        fit_prms([1,3,5]) = [mean(G_o,'omitnan'),m,a];
end

% Grab initial parameter set and bounds if needed
x_i = x_i(prmidx);
lb = Fit_opts.x_lb(prmidx);
ub = Fit_opts.x_ub(prmidx);

% Set weights for inversion
switch Fit_opts.WeightFlag
    case 0 % Even weights
        W = ones(size(In,1),1);
    case 1 % Inverse weighting
        W = 1./In.(Fit_opts.WeightVar);
        W(isnan(W)) = 0;
    case 2 % Inverse squared weighting
        W = 1./In.(Fit_opts.WeightVar).^2;
        W(isnan(W)) = 0;
end

% 2. FIT MODEL
%'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

switch Fit_opts.SolverFlag
    case 1 % fitnlm, nonlinear least squares
        
        % Setup fit options
        opts = statset('nlinfit');
        if Fit_opts.Robust == 1
            opts.Robust = 'on';
            opts.RobustWgtFun = Fit_opts.RobustWgtFun;
        else
        end
        
        % Fit
        fit = fitnlm(In,mdl,x_i,'PredictorVars',predvars,...
            'ResponseVar',respvar,'Options',opts);
        
        % Store outputs
        x_opt = fit.Coefficients.Estimate;
        fit_est = fit.Fitted;
        AIC = fit.ModelCriterion.AIC;
        BIC = fit.ModelCriterion.BIC;
        covmat = fit.CoefficientCovariance;
        
    case {2,3,4} % Use older solvers that need more adjustments
        
        % Create array of input data for older solvers
        In_array = table2array(In(:,predvars));
        switch Fit_opts.SolverFlag
            case 2 % lsqcurvefit, constrained NLS (no weighting allowed)
                
                % Setup fit options
                opts = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
                opts.Display = 'none';
                
                % Fit
                [x_opt,~,res,~,~,~,J]  = lsqcurvefit(mdl,x_i,In_array,...
                    In.(respvar),lb,ub,opts);
                
                % Calculate some outputs
                covmat = full(inv(J'*J))*...
                    (res'*res)./(length(res) - length(x_opt));
                fit_est = In.(respvar) + res;
                AIC = NaN;
                BIC = NaN;
                
            case {3,4} % More general solver options
                
                % Specify loss function
                if contains(Fit_opts.LossFxnFlag,'L1') % L1 norm
                    lossfxn = @(x) sum(abs(In.(respvar) - mdl(x,In_array)).*W);
                elseif contains(Fit_opts.LossFxnFlag,'L2') % L2 norm
                    lossfxn = @(x) sum((In.(respvar) - mdl(x,In_array)).^2.*W);
                elseif contains(Fit_opts.LossFxnFlag,'LCE') % Least cost
                    lossfxn = @(x) 1 - calcLCE(In.(respvar),mdl(x,In_array));
                end
                
                switch Fit_opts.SolverFlag
                    case 3 % fminsearch, Nelder-Mead simplex method
                        
                        % Update the the loss fxn
                        opts = optimset('fminsearch');
                        opts.Display = 'none';
                        x_opt  = fminsearch(lossfxn,x_i,opts);
                        fit_est = mdl(x_opt,In_array);
                        covmat = zeros(length(x_opt));
                        AIC = NaN;
                        BIC = NaN;
                        
                    case 4 % fmincon,constrained optimization
                        
                        % Update the the loss fxn
                        opts = optimoptions('fmincon','Algorithm','interior-point');
                        opts.Display = 'none';
                        [x_opt,~,~,~,~,hessian]  = fmincon(lossfxn,x_i,[],[],[],[],lb,ub,[],opts);
                        fit_est = mdl(x_opt,In_array);
                        covmat = 2*var(In.(respvar) - fit_est)/hessian;
                        AIC = NaN;
                        BIC = NaN;
                        
                end
        end
end

% Store remainder of outputs
fit_prms(prmidx) = x_opt;
fit_mdl.fit_mdl = @(x) mdl(x_opt,x);
fit_mdl.NumCoefficients = length(x_opt);
fit_mdl.predvars = predvars;
fit_mdl.covmat = covmat;
fit_mdl.AIC = AIC;
fit_mdl.BIC = BIC;

% Extra step if VPDFlag = 2
if isequal(Fit_opts.VPDFlag,2)
     In_array = table2array(In(:,predvars));
    [~,~,~,VPD_l] = calcGcPMforFit(fit_prms(1),fit_prms(2),fit_prms(3),...
        fit_prms(4:end),In_array,Fit_opts);
else
    VPD_l = [];
end

end
