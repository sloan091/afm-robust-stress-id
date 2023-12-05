function [Obskeep,GPPkeep,GPPdays] = filtGPP(PP_opts,Obs)

        % Grab 95th perctile of daily GPP data
        tmp = load([PP_opts.Site_ID,'_Daily']);
        D = tmp.Obs_Daily;
        GPPmax = prctile(D.GPP_NT_VUT_REF,95);
        
        % Smooth with 15 average filter
        GPPsm = smoothdata(D.GPP_NT_VUT_REF,'movmean',days(15),'omitnan',...
            'SamplePoints',D.Time);
        
        % Apply threhsold
        switch PP_opts.GPPFiltFlag
            case 1 %  > 10% of 95th GPP prctile
                GPPdays = GPPsm > 0.1*GPPmax;
            case 2 % > 50% of 95th GPP prctile
                GPPdays = GPPsm > 0.5*GPPmax;
        end
        
        % Remove half-hourly observations not on those days
        GPPdays = D.Time(GPPdays);
        hhtime = dateshift(Obs.Time,'start','day');
        GPPkeep = ismember(hhtime,GPPdays); 
        Obskeep = Obs(GPPkeep,:);

end