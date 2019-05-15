
% --------------------------------------------------------------------
% visualisation subroutine
% --------------------------------------------------------------------


% clean up
clc;
clear all;
close all;

Ncolors  = 25;


filename = 'Stats_G_trend_noname.xls';
T = xlsread(filename);

lon = T(:,1);
lat = T(:,2);
Gmean = T(:,9);

Gmin  = -2;
Gmax  = 2;
col   = floor( min( max( (Gmean-Gmin)/(Gmax-Gmin)-1e-6 ,0) ,1)*(Ncolors-1)+1);


cmap = flipud(cbrewer('div','RdBu',Ncolors));

 
figure
set(gcf, 'Position', [100, 100, 2000, 1000])
m_proj('robinson','lat',[-60 90]); 
% m_proj('miller','lat',82);            
m_coast('color',[0.1 0.1 0.1]); 
% m_grid('linestyle','none','box','fancy','tickdir','out');
m_grid('linestyle','none','tickdir','out','linewidth',1);
indices = find(~isnan(Gmean));
for i=1:length(indices)
    ind = indices(i);
%     if(ind==95)
%        'test' 
%     end
    m_line(lon(ind), lat(ind), 'linestyle', '-', 'marker', 'o', 'markersize', 9, 'MarkerFaceColor', cmap(col(ind),:), 'color', [0.8 0.8 0.8]);
end
colormap(cmap);
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0, 1, Ncolors/2) ; %Create 8 ticks from zero to 1
cbh.TickLabels = floor(100*linspace(Gmin, Gmax, Ncolors/2))/100 ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
title('Mean G FLUXNET (W/m^2)','Fontsize',24)
set(gca,'fontsize',24) 
% export_fig pct_irr_lowres -transparent -png % save figure
export_fig Figures/mean_G -m10 -transparent -png% save figure high quality

'end script'