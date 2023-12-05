% =========================================================================
% Name   : importFluxnetData.m
% Author : Brandon Sloan
% Date   : 6/14/21
%
% DESCRIPTION
% This function reads in FLUXNET2015 data for a select site from its native
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
function [Obs,pname] = importFluxnetDataDaily(folder,tres)

% Grab file name
fname = dir([folder,'\*_FULLSET_',tres,'_*']);
path = [folder,'\',fname.name];
pname = fname.name(5:10);

% Import file with selected variables
opts = detectImportOptions(path);
opts = setvaropts(opts,'TreatAsMissing','-9999');
tempv = opts.VariableNames;
rem1 = contains(tempv,{'_JSB','_ERA','_CUT'});
rem2 = contains(tempv,{'RECO_','GPP_','NEE_'}) & ~contains(tempv,'_REF');
keep  = ~rem1 &  ~rem2;
opts.SelectedVariableNames = tempv(keep);
dt = readtable(path,opts);

% Convert to MATLAB timetable
dates = num2str(dt.TIMESTAMP);
dt.Time = datetime(datevec(dates,'yyyymmdd'));
Obs = table2timetable(dt,'RowTimes','Time');

% Update name
pname(3)='_';

end

