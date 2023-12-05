function LCE = calcLCE(Obs,Model)

% GOF columns: RMSE,NSE,Pbias,R,cRMSE,StdRatio,MeanRatio,MeanET,AIC,BIC
cor = corrcoef(Obs,Model,'Rows','complete');
r = cor(2);
a = std(Model,'omitnan')/std(Obs,'omitnan');
b = mean(Model,'omitnan')/mean(Obs,'omitnan');
LCE = 1 - sqrt((r.*a - 1).^2 + (r./a - 1).^2 + (b - 1).^2);

end

