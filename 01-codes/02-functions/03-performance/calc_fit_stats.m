function  [RMSE,NSE,Pbias,R,cRMSE,StdRatio,MeanRatio,nRMSE,nNSE] = calc_fit_stats(Obs,Mdl)

RMSE = sqrt(nanmean((Mdl - Obs).^2));
nRMSE = RMSE./nanmean(Obs);
NSE = 1 - nansum((Mdl - Obs).^2)./nansum((Obs - nanmean(Obs)).^2);
nNSE = 1./(2 - NSE);
Pbias = nansum(Mdl - Obs)/nansum(Obs);
Rtmp = corrcoef(Mdl,Obs,'Rows','pairwise');
R = Rtmp(2);

cRMSE = sqrt(nanmean(((Obs - nanmean(Obs))-(Mdl - nanmean(Mdl))).^2));
StdRatio = nanstd(Mdl)/nanstd(Obs);
MeanRatio = nanmean(Mdl)/nanmean(Obs);

end