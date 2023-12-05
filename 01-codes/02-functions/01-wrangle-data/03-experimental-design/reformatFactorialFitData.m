% This function reformats my fit data from 3D matrices to nice tables with
% headings.
function [prms_trn,gof_trn,gof_tst,info_trn,info_tst] = ...
    reformatFactorialFitData(prms_trn,gof_trn,gof_tst,info_trn,info_tst,binFlag,nbins,nboot)

switch binFlag
    case 1 % Binned fits
        % Reformat output files to tables
        prmvn = {'ExpID','BinID','G_o','G_1','m',...
            'G_1_VPD_m','G_v'};
        gofvn = {'ExpID','BinID','RMSE','NSE','Pb','R','cRMSE','StdRatio',...
            'MeanRatio','MeanFlux','AIC','BIC'};
        infovn = {'ExpID','BinID','InfoID','Ap','AfT','Afu1','Afu2','Afs','Afr','Afp'};
        BinID = reshape(repmat((1:nbins),...
            nboot,size(prms_trn,3)),[],1);
        BinIDinfo = repmat(reshape(repmat(1:size(info_trn,3),...
            size(info_trn,1),1),[],1),size(info_trn,4),1);
        ExpID = reshape(repmat(1:size(prms_trn,3),...
            size(prms_trn,1),1),[],1);
        ExpIDinfo = reshape(repmat(1:size(info_trn,4),size(info_trn,1).*...
            size(info_trn,3),1),[],1);
        hlpr = @(x,vn) array2table([ExpID,BinID,conv3Dto2D(x)],'VariableNames',vn);
        hlprinfo = @(x,vn) array2table([ExpIDinfo,BinIDinfo,conv4Dto2D(x)],'VariableNames',vn);
        prms_trn = hlpr(prms_trn,prmvn);
        gof_trn = hlpr(gof_trn,gofvn);
        gof_tst = hlpr(gof_tst,gofvn);
        info_trn = hlprinfo(info_trn,infovn);
        info_tst = hlprinfo(info_tst,infovn);
        
        % Make bins categorical
        cathelp = @(T) convertvars(T,{'ExpID','BinID'},'categorical');
        prms_trn= cathelp(prms_trn);
        gof_trn = cathelp(gof_trn);
        gof_tst = cathelp(gof_tst);
        info_trn = cathelp(info_trn);
        info_tst = cathelp(info_tst);
        
    case 0
%         % Reformat output files to tables
%         prmvn = {'G_o','G_1','m',...
%             'G_1_VPD_m','G_v'};
%         gofvn = {'RMSE','NSE','Pb','R','cRMSE','StdRatio',...
%             'MeanRatio','MeanFlux','AIC','BIC'};
%         infovn = {'InfoID','Ap','AfT','Afu1','Afu2','Afs','Afr','Afp'};
%         hlpr = @(x,vn) array2table(x,'VariableNames',vn);
%         prms_trn = hlpr(prms_trn,prmvn);
%         gof_trn = hlpr(gof_trn,gofvn);
%         gof_tst = hlpr(gof_tst,gofvn);
%         info_trn = hlpr(info_trn,infovn);
%         info_tst = hlpr(info_tst,infovn);
%         info_trn.InfoID = categorical(info_trn.InfoID);
%         info_tst.InfoID = categorical(info_tst.InfoID);
        
end
% Remove AIC and BIC from gof_tst
gof_tst = removevars(gof_tst,{'AIC','BIC'});

end