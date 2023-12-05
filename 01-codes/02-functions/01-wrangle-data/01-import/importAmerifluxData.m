% =========================================================================
% Name   : importAmerifluxData.m
% Author : Brandon Sloan
% Date   : 6/14/21
%
% DESCRIPTION
% This function reads in Ameriflux data for a select site from its native
% .csv format and creates a matlab timetable.
%
% INPUTS
%   folder - The folder path of FLUXNET2015data
%   tres   - Indicator for temporal resoultion
%
% OUTPUTS
%   Obs - MATLAB timetable of observation
%
% =========================================================================
function Obs = importAmerifluxData(folder,tres)

% Begin function
fname = dir([folder,'*_',tres,'_*']);
path = [folder,fname.name];
pname = fname.name(1:10);

% These variables are not relevant to our current study and are removed
vars2del = {''};
% vars2del = {'SW_IN_POT','SW_IN_F','SW_IN_F_QC','LW_IN_F','LW_IN_F_QC','WD'...
%     'PPFD_IN','CO2_F_MDS','CO2_F_MDS_QC','LE_CORR_25','LE_CORR_75','H_CORR_25','H_CORR_75'...
%     'NEE_VUT_REF','NEE_VUT_REF_QC','NEE_VUT_REF_RANDUNC','NEE_VUT_25','NEE_VUT_50','NEE_VUT_75'...
%     'NEE_VUT_25_QC','NEE_VUT_50_QC','NEE_VUT_75_QC','RECO_NT_VUT_REF','RECO_NT_VUT_25','RECO_NT_VUT_50'...
%     'RECO_NT_VUT_75','GPP_NT_VUT_REF','GPP_NT_VUT_25','GPP_NT_VUT_50','GPP_NT_VUT_75','RECO_DT_VUT_REF'...
%     'RECO_DT_VUT_25','RECO_DT_VUT_50','RECO_DT_VUT_75','GPP_DT_VUT_REF','GPP_DT_VUT_25','GPP_DT_VUT_50','GPP_DT_VUT_75','RECO_SR','RECO_SR_N'};
opts = detectImportOptions(path);
opts = setvaropts(opts,'TreatAsMissing','-9999');
tempv = opts.VariableNames;
Lia = ismember(tempv,vars2del);
opts.SelectedVariableNames = tempv(~Lia);
dt = readtable(path,opts);

if strmatch(tres','HH')
    
    % This portion converts everything to MATLAB datenumbers
    dates = num2str(dt.TIMESTAMP_START);
    dt.TIMESTAMP_START = datetime(datevec(dates,'yyyymmddHHMM'));
    dates = num2str(dt.TIMESTAMP_END);
    dt.TIMESTAMP_END = datetime(datevec(dates,'yyyymmddHHMM'));
    
else
    
    % This portion converts everything to MATLAB datenumbers
    dates = num2str(dt.TIMESTAMP);
    dt.TIMESTAMP_START = datetime(datevec(dates,'yyyymmdd'));

end

% Need to replace hyphen in case I want to store as a structure with this
% name
pname(7)='_';

% Convert to MATLAB timetable
Obs = table2timetable(dt);

end

