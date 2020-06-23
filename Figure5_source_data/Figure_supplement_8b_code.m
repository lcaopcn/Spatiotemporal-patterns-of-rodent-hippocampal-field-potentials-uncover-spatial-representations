clear; clc;
warning('off');

%% Load and PreProcess Data
err = load('Figure_supplement_8b_results_data.mat');

figure('Color','white');

clo = {[1 0 0],[0 0 1],[0 1 0]};
subplot(1,2,2)
h(:,1) = cdfplot(err.err_nd_All(:,1)); hold on;
h(:,2) = cdfplot(err.err_nd_All(:,2)); xlim([0,50]); title('');
grid off; yticks([0,0.5,1]); xticks([0,25,50]);
xlabel(''); ylabel('');
xlabel('Decoding error (cm)'); ylabel('CDF');

set( h(:,1), 'LineStyle', ':', 'LineWidth',1.5, 'Color', clo{3});
set( h(:,2), 'LineStyle', '-', 'LineWidth',1.5, 'Color', clo{3});   

trackLen = err.trackLen;
edges = linspace(0,2*pi,trackLen/2);
e{1}=edges; e{2}=edges;
ex = edges+mean(diff(edges))/2;
[mx,my] = meshgrid(ex,ex);
posbinc = 1:2:trackLen-1;

cmaxrange = 0.1;

subplot(1,2,1)
colormap(flipud(bone));
[n,~]=hist3([err.rpos*2*pi,err.pos_hat(:,2)],'Edges',e);
ns = imgaussfilt(n/4,2);
imagesc(posbinc,posbinc,ns); set(gca,'YDir','normal');
axis image; clim([0,cmaxrange]);
xlim([0 1]*trackLen);
ylim([0 1]*trackLen); 
xlabel('True position (cm)'); ylabel('Estimated (cm)');
set(gca,'TickDir','out')
set(gca,'XTick',[0 trackLen])
set(gca,'YTick',[0 trackLen])
ax = gca; colormap(ax,flipud(bone)); clim([0,cmaxrange]);
cb = colorbar(ax);
cb.Label.String = "Density (count/cm^2)"; cb.Label.Position = [1,cmaxrange/2];
set(cb,'XTick',[0,cmaxrange]); cb.FontSize = 8;

savefig(gcf,'Figure_supplement_8b','compact');
