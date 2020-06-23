clear; clc;
warning('off');

%% Load and PreProcess Data
load('Figure2_results_data_decoding');

figure('Color','white');
edges = linspace(0,2*pi,trackLen/2);
e{1}=edges; e{2}=edges;
ex = edges+mean(diff(edges))/2;
[mx,my] = meshgrid(ex,ex);
posbinc = 1:2:289;

cmaxrange = 0.5;
for ib=1:2
    subplot(2,1,ib)
    colormap(flipud(bone));
    [n,~]=hist3([rpos*2*pi,pos_hat(:,ib)],'Edges',e);
    ns = imgaussfilt(n/4,1);
    imagesc(posbinc,posbinc,ns); set(gca,'YDir','normal');
    axis image; clim([0,cmaxrange]);
    xlim([0 1]*trackLen); 
    ylim([0 1]*trackLen); 
    if ib == 2, xlabel({'Position from      (cm)'}); end
    ylabel(['Position from      (cm) '])
    set(gca,'TickDir','out')
    set(gca,'XTick',[0 trackLen])
    set(gca,'YTick',[0 trackLen])
end
ax = axes; colormap(ax,flipud(bone)); clim([0,cmaxrange]);
cb = colorbar(ax); cb.Position = [0.85,0.25,0.05,0.5]; 
cb.Label.String = "Density (count/cm^2)"; cb.Label.Position = [1,cmaxrange/2];
set(cb,'XTick',[0,cmaxrange]); cb.FontSize = 8; ax.Visible = 'off';

savefig(gcf,'Figure2c','compact');
