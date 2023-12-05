function  [Hstruct,Istruct,Partstruct] = calcInfoDiagnostics(x1raw,x2raw,yraw,nbins)

% Rescale data from 0 to 1
x1 = rescaleData(x1raw);
x2 = rescaleData(x2raw);
y  = rescaleData(yraw);

% Create equidistant bins
bins = repmat(linspace(0,1,nbins)',1,3);

% Calculate all needed information entropies
[Hx1,Hx2,Hy,Hx1x2,Hx1y,Hx2y,Hx1x2y] = calc_Hx_combinations(x1,x2,y,bins);
Hstruct = v2struct(Hx1,Hx2,Hy,Hx1x2,Hx1y,Hx2y,Hx1x2y);

% Calculate scaled redundant information
Ix1_x2 = Hx1 + Hx2 - Hx1x2;
Ix1_y  = Hx1 + Hy - Hx1y;
Ix2_y  = Hx2 + Hy - Hx2y;
Ix1_ygx2 = Hx1x2 + Hx2y - Hx1x2y - Hx2;
Ix1_x2_y = Ix1_ygx2 - Ix1_y;
Rmin = max([0,-Ix1_x2_y]);
Rmmi = min([Ix1_y,Ix2_y]);
scaledR = Ix1_x2/min([Hx1,Hx2]);
Rs = scaledR*(Rmmi - Rmin) + Rmin;


% Calculate synergistic information
S = Ix1_x2_y + Rs;

% Calculate unique information
Ux1 = Ix1_y + Rs;
Ux2 = Ix2_y + Rs;

% Gather components
partition  = [Ux1,Ux2,S,Rs];

% Check partition is correct
Ix1x2_y = Hx1x2 + Hy - Hx1x2y;
Residual = (Ix1x2_y - sum(partition))/Ix1x2_y;
Istruct = v2struct(Ix1_x2,Ix1_y,Ix2_y,Ix1_ygx2,Ix1_x2_y,Ix1x2_y,Residual);

% Normalize parition
nrmpart = partition/sum(partition);
% nrmpart = partition/Ix1x2_y;

Partstruct = v2struct(Rs,S,Ux1,Ux2,partition,nrmpart);

end


function xrs = rescaleData(x1)

xrs = (x1 - min(x1))./(max(x1) - min(x1));

end

function px = calc_px(x,bins)

d = size(bins,2);
n = size(x,1);

switch d
    case 1
        hct = histcounts(x,bins);
        px = hct/n;
    case 2
        hct = histcounts2(x(:,1),x(:,2),bins(:,1),bins(:,2));
        px = hct/n;
    case 3
        [hct,a,b,c] = histcn(x,bins(:,1),bins(:,2),bins(:,3));
        px = hct/n;
end

end

function [px1,px2,py,px1x2,px1y,px2y,px1x2y] = calc_px_combinations(x1,x2,y,bins)

% Discrete probabilities
px1 = calc_px(x1,bins(:,1));
px2 = calc_px(x2,bins(:,2));
py  = calc_px(y,bins(:,3));

% Discrete 2D joint probabilities
px1x2 = calc_px([x1,x2],[bins(:,1),bins(:,2)]);
px1y  = calc_px([x1,y],[bins(:,1),bins(:,3)]);
px2y  = calc_px([x2,y],[bins(:,2),bins(:,3)]);

% Discrete 3D joint probabilities
px1x2y = calc_px([x1,x2,y],bins);


end

function [Hx1,Hx2,Hy,Hx1x2,Hx1y,Hx2y,Hx1x2y] = calc_Hx_combinations(x1,x2,y,bins)

[px1,px2,py,px1x2,px1y,px2y,px1x2y] = calc_px_combinations(x1,x2,y,bins);

% Information Entropy
Hx1 = calcEntropy(px1);
Hx2 = calcEntropy(px2);
Hy  = calcEntropy(py);

% Joint 2D Information Entropy
Hx1x2 = calcEntropy(px1x2);
Hx1y  = calcEntropy(px1y);
Hx2y  = calcEntropy(px2y);

% Joint 3D Information Entropy
Hx1x2y = calcEntropy(px1x2y);

end

function Hx = calcEntropy(px)

ip = px.*log2(px);
Hx = -nansum(nansum(ip(:)));

end