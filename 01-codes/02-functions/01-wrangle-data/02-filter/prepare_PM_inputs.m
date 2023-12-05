% =========================================================================
% Name   : prepare_PM_inputs.m
% Author : Brandon Sloan
% Date   : 2/22/22
%
% DESCRIPTION
%
% =========================================================================

function PM_inp = prepare_PM_inputs(PM_inp,Obs)

% Unit corrections for PM
Obs.RH = Obs.RH/100;
Obs.P_a = Obs.P_a*1000;

% Create PM input dataset
PM_inp = PM_inp;
PMnames = Obs.Properties.VariableNames;
for ii = 1:length(PMnames)
    PM_inp.(PMnames{ii}) = Obs.(PMnames{ii});  
end
PM_inp.time = Obs.Time;
PM_inp.G(isnan(PM_inp.G)) = 0;

end