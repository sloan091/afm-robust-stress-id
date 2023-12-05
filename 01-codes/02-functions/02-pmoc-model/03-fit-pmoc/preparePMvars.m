function Out = preparePMvars(In,PM_opts)
    
    % Create table Penman_Monteith code
    G_co = In.G_c;
    In = timetable2table(In);
    In = In(:,2:end);
    
    % Correct FLUXNET observation units for Penman-Monteith
    In.RH = In.RH/100;
    In.P_a = In.P_a*1000;
    In.G(isnan(In.G)) = 0;
    In.G_c = ones(size(In,1),1);
    
    % Estimate Penman-Monteith quantities using G_c = 1
    [~,~,~,PM] = calc_ET_Penman_Monteith(In,PM_opts);
    
    % Reset G_c to observed values
    In.G_c = G_co;
    
    % Final variables
    Out = [In,PM];

end