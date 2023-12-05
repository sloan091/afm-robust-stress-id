% =========================================================================
% Name   : step_02_export_fit_to_excel_for_r.m
% Author : Brandon Sloan
% Date   : 5/31/22
%
% DESCRIPTION
% This filters and reformats the PMOC model fits and exports them to .xlsx
% files that can be read into R for the robustness analysis in Sloan and
% Feng (2023).
%   
% =========================================================================

clc
clear
close all
load pmoc_default_opts.mat
load final_ec_site_properties.mat

% Print name
basefp = '.\03-outputs\01-estimated-pmoc-parameters\01-matlab\';
svfp = '.\03-outputs\01-estimated-pmoc-parameters\02-excel\';

% Variables
vn = {'SWC_nep','SubTrtID','G_o','G_1','m','G_1_VPD_m','G_v','T_ET',...
    'nRMSE','nNSE','Pb','R','cRMSE','StdRatio','MeanRatio','MeanFlux',...
    'LCE','Ap','AfT','Afp'};
keep = {'SWC_nep','SubTrtID','G_o','G_1','m','G_1_VPD_m','G_v','T_ET',...
    'nRMSE','LCE','Ap','AfT','Afp'};

% Absolute outlier threshold for G1 and G1/VPD^m
outThresh = [-10,50];

% Select sites
files = extractfield(dir([basefp,'\*.mat']),'name');
nsites = length(files);

% Create experimental design matrix
ED = fullfact([2,2,4,3,4,3,2,2]);
fctnames =  append('Fct',strsplit(num2str(1:size(ED,2))));
ED = array2table(ED,'VariableNames',fctnames);
ED.ExpID = categorical((1:size(ED,1)).');

% Catch badfiles
badfiles = {};
remPrc = NaN(nsites,2);

% Cycle through sites
for jj = 1:nsites
    
    try
        % Set site properties
        filejj = files{jj};
        sitenmjj = filejj(1:6);
        
        % Load results
        load([basefp,filejj]);
        
        % Update GOF and info and add to parameters
        gof_tst.nRMSE = gof_tst.RMSE./gof_tst.MeanFlux;
        gof_tst.LCE =  1 - sqrt((gof_tst.R.*gof_tst.StdRatio - 1).^2 +...
            (gof_tst.R./gof_tst.StdRatio  - 1).^2 + ...
            (gof_tst.MeanRatio  - 1).^2);
        info_keep = info_tst(info_tst.InfoID == 1,'Ap');
        info_keep.AfT1 = info_tst{info_tst.InfoID == 1,'AfT'}; 
        info_keep.AfT2 = info_tst{info_tst.InfoID == 2,'AfT'}; 
        info_keep.AfT3 = info_tst{info_tst.InfoID == 3,'AfT'};
        R = [prms_trn,gof_tst(:,{'NSE','LCE'}),info_keep];
        
        % Add SWC_nep
        bintab = table(categorical(1:length(Fit_opts.binedges) - 1)',...
            movmean(Fit_opts.binedges,2,'Endpoints','discard')',...
            'VariableNames',{'BinID','SWC_nep'});
        R = outerjoin(R,bintab,'Type','Left','Keys','BinID','MergeKeys',true);
        R = outerjoin(R,ED,'Type','Left','Keys','ExpID','MergeKeys',true);
        
        % Remove outliers and Nan
        [R_Table,outPrc,nanPrc] = OutlierRemoval(R,'BinID',outThresh);
        remPrc(jj,:) = [outPrc,nanPrc];
        
        % Remove and re-arrange variables and save
        R_Table = movevars(R_Table,[fctnames,'SWC_nep'],...
            'Before','ExpID');
        R_Table = removevars(R_Table,{'BinID'});
        writetable(R_Table,[svfp,sitenmjj,'_FF_AFM','.csv']);
        
    catch ME
        badfiles = [badfiles;[sitenmjj,': ',ME.message]];
        
    end
end

% Save the outlier percentage for each site
remPrc = array2table(remPrc*100, "VariableNames",{'Outliers', 'NaN'});
remPrc.Site_ID = SiteProp.Site_ID;
save('.\01-codes\01-scripts\02-fit-pmoc-model\factorial_fit_outliers','remPrc')

function [Out,outPrc,nanPrc] = OutlierRemoval(In,binName,outThresh)

binID = unique(In.(binName));

% Check which SWC bins must be excluded due to inadequate treatment level
% coverage
FctNames = ["Gc","Flx","Prm","VPD","Fit","SEB","GPP","LAI"];
chkmat = [];
FctLvl = [];
for ii = 1:8
    fctii = ['Fct',num2str(ii)];
    nanct_ii = grpstats(In(:,{'BinID',fctii,'G_v'}),{'BinID',fctii},...
        @(x) sum(~isnan(x)),'VarNames',{'BinID','Fct','ignore','chk'});
    lvl = unique(nanct_ii.Fct)';
    FctLvl = [FctLvl,FctNames(ii) + lvl];
    chkmat = [chkmat,reshape(nanct_ii.chk,length(lvl),10)'];
end

% Remove bins that were not run
allFmiss = sum(chkmat,2) == 0;
In = In(ismember(In.BinID,binID(~allFmiss)),:);

% Remove very large values
keep1 = ~isnan(In.G_1);
keep2 = In.G_1 > outThresh(1) & In.G_1 < outThresh(2) & ...
    In.G_1_VPD_m > outThresh(1) & In.G_1_VPD_m < outThresh(2) & ...
    In.LCE >= -10;
Out = In(keep1 & keep2,:);

% Re-sort the data
Out = sortrows(Out,{'ExpID',binName});

% Store removal percentage
nanPrc = sum(~keep1)/size(In,1);
outPrc = sum(~keep2 & keep1)/size(In,1);

end


