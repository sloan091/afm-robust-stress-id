% =========================================================================
% Name   : myBootstrapwTest.m
% Author : Brandon Sloan
% Date   : 7/13/22
%
% DESCRIPTION
% This function performs bootstraps parameter fits (prm) and goodness of
% fit (gof) statistics by training (_trn) a model and testing (_tst) it 
% against out-of-sample observations.
% =========================================================================

function [prms_trn,gof_trn,gof_tst,info_trn,info_tst] = ...
    train_fit(nboot,trainfxn,testfxn,inputs,ParFlag,nprms,nstats)

% Pre-allocate
prms_trn = NaN(nboot,nprms);
gof_trn = NaN(nboot,nstats);
gof_tst = gof_trn;
info_trn = NaN(nboot,8*3);
info_tst = info_trn;
nall = size(inputs,1);

switch ParFlag % Keep ParFlage at 0 since parallelization happens earlier.
    case 0 % Perform bootstrap sequentially
        
        for ii = 1:nboot
            
            % Create bootstrap random sample
            %rng(0,'twister');
            [trn_set,trnidx] = datasample(inputs,nall,1,'Replace',true);
            tstidx = find(~ismember(1:nall,unique(trnidx)));
            
            % Split bootstrap sample into training and test set
            tst_set = inputs(tstidx,:);
            
            try
                % Train model to data
                [prms_trn(ii,:),gof_trn(ii,:),trn_mdl,info_trn(ii,:)] =...
                    trainfxn(trn_set);
                
                % Test the trained model on out-of-sample observations
                [gof_tst(ii,:),info_tst(ii,:)] = testfxn(trn_mdl,tst_set);

            catch
            end
            
        end
    case 1 % Perform bootstrap in parallel
        parfor ii = 1:nboot
            
            % Create bootstrap random sample
            [trn_set,trnidx] = datasample(inputs,nall,1,'Replace',true);
            tstidx = find(~ismember(1:nall,unique(trnidx)));
            
            % Split bootstrap sample into training and test set
            tst_set = inputs(tstidx,:);
            
            try
                % Run training and test functions
                [prms_trn(ii,:),gof_trn(ii,:),trn_mdl,info_trn(ii,:)] =...
                    trainfxn(trn_set);
                
                % Test the trained model on out-of-sample observations
                 [gof_tst(ii,:),info_tst(ii,:)] = testfxn(trn_mdl,tst_set);
            catch
            end
            
        end
end

% Reshape info matrix
info_trn = sortrows(reshape(info_trn',8,[])',1);
info_tst = sortrows(reshape(info_tst',8,[])',1);

     
end