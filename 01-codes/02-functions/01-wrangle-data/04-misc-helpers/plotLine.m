function [h,gax] = plotLine(xdata,ydata,xlab,ylab,DN,LS,LW,...
    axLW,xlim,ylim,fontsize,fontname,color,position,varargin)

if isempty(varargin)
    gax = gca;
else
    gax = varargin{:};
end

set(groot,'defaulttextinterpreter','tex');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(groot, 'defaultLegendInterpreter','tex');
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h = plot(xdata,ydata,'Color',color,'DisplayName',DN,'LineStyle',LS,...
    'Linewidth',LW,'Marker','none','Parent',gax);

if isempty(position)
set(gax,'LineWidth',axLW,'ycolor','k','xcolor','k')    
else
set(gax,'LineWidth',axLW,'ycolor','k','xcolor','k',...
    'ActivePositionProperty','position','units','centimeters','Position',position);
end

if ~isempty(xlim)
    set(gax,'xlim',xlim);
else
end

if ~isempty(ylim)
    set(gax,'ylim',ylim);
else
end

xlabel(gax,xlab);
ylabel(gax,ylab);
set(gax,'TickDir','out','FontName',fontname,'FontUnits','point',...
    'FontSize',fontsize,...
    'LabelFontSizeMultiplier',1)
box(gax,'off');
f = gax.Parent;
if strmatch(class(f),'matlab.graphics.layout.TiledChartLayout')
    f = f.Parent;
else
end

set(f,'Color','w');
if isempty(position)
    
else
f.Units = 'centimeters';
f.Position(3) = position(3) + 2*gax.TightInset(1);
f.Position(4) = position(4) + 2*gax.TightInset(2);
end

end

