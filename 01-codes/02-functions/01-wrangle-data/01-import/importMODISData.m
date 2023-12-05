% =========================================================================
% Name   : importMODISData.m
% Author : Brandon Sloan
% Date   : 10/9/21
%
% DESCRIPTION
% This function reads in MODIS LAI/FPAR for a select site from its native
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
function Obs = importMODISData(folder,type)

% Begin function

fnames = extractfield(dir([folder]),'name');

if contains(type,{'NDVI','EVI'})
    fname = regexp(fnames, regexptranslate('wildcard','stat*13*'),'match','once');
    tst = @(x) extractAfter(x,find(x == '_',1,'last'));
elseif contains(type,{'LAI','PAR'})
    fname = regexp(fnames, regexptranslate('wildcard','stat*15*'),'match','once');
    tst = @(x) extractBefore(x,find(x == '_',1,'first'));
else
    error('Invalid File Type');
end

fname(cellfun('isempty',fname)) = [];

% Create 3 separate paths
% path1 = [folder,'\',fname(4).name];
% path2 = [folder,'\',fname(3).name];
path = [folder,'\',fname{:}];

% Grab scaled and QC applied fPAR/LAI
opts = detectImportOptions(path);
opts = setvaropts(opts,'TreatAsMissing','NaN');
dt = readtable(path,opts);

time = datetime(dt.dt,'InputFormat','dd/MM/uuuu');
vns = unique(dt.band);
vn = cellfun(tst,vns,'UniformOutput',0);
% for i = 1:length(time)
%     year = str2num(time{i,1}(2:5));
%     doy = str2num(time{i,1}(6:end));
%     dates(i,1) = datetime(year,1,doy);
% end

% fPAR data
v1_site = dt.value_center(1:2:end,1);
v1_mean = dt.value_mean(1:2:end,1);
v1_std = dt.value_standard_deviation(1:2:end,1);

% LAI data
v2_site = dt.value_center(2:2:end,1);
v2_mean = dt.value_mean(2:2:end,1);
v2_std = dt.value_standard_deviation(2:2:end,1);

% Final dataset
Obs = timetable(unique(time),v1_site,v1_mean,v1_std,...
    v2_site,v2_mean,v2_std);
lbl = upper(string(vn)) + string({'_site','_mean','_std'});
Obs.Properties.VariableNames = lbl';

end

