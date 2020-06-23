clear; clc
warning('off');

load('.\Figure_supplement_6_source_data.mat');
[W,H] = nnmf(rFPAuhf',64);
H = zscore(H,0,1);
rH = int64(resample2poission(H',10,rand()));

makers = {'s','^','o','h'};
cdfval = [0.1,0.3,0.6,0.9];

%% plot
fig = figure('Color','w');
subplot(2,2,1)
h = histogram(H); h.Normalization = 'pdf';
h.NumBins = 100; h.FaceColor = 'b'; h.FaceAlpha = 0.7;
xlim([-2,4]); h.EdgeColor = 'w'; h.EdgeAlpha = 0.7; 
ylabel('PDF','FontSize',12); xlb = xlabel('$y$','FontSize',12); set(xlb,'Interpreter','latex');
title('Original distribution');

subplot(2,2,3)
Hl = reshape(H,[],1);
[f,x] = ecdf(Hl); plot(x,f,'LineWidth',2);
hold on; box on;
ylabel('CDF','FontSize',12); xlim([-2,4]);
for isc = 1:length(cdfval)
    idx = find(f(1:end-1)<cdfval(isc) & f(2:end)>cdfval(isc));
    scatter(x(idx),f(idx),72,'r',makers{isc},'filled');
end
xlb = xlabel('$y$','FontSize',12); set(xlb,'Interpreter','latex');

subplot(2,2,2)
h = histogram(rH); h.Normalization = 'pdf';
h.BinWidth = 1; h.FaceColor = 'b'; h.FaceAlpha = 0.7;
h.EdgeColor = 'w'; h.EdgeAlpha = 0.7; 
xlim([0,23]); ylabel('PDF','FontSize',12); xlb = xlabel('$\hat y$','FontSize',12); set(xlb,'Interpreter','latex');
title('Target distribution (Poissanoon)');

subplot(2,2,4)
rHl = reshape(rH,[],1);
[f,x] = ecdf(double(rHl)); plot(x,f,'LineWidth',2);
hold on; box on;
ylabel('CDF','FontSize',12); xlim([0,23])
for isc = 1:length(cdfval)
    idx = find(f(1:end-1)<cdfval(isc) & f(2:end)>cdfval(isc))+1;
    scatter(x(idx),f(idx),72,'r',makers{isc},'filled');
end
xlb = xlabel('$\hat y$','FontSize',12); set(xlb,'Interpreter','latex');

annotation('arrow', [.4 .4], [.65 .4],'Linewidth',1.5);
annotation('arrow', [.65 .65], [.65 .4],'Linewidth',1.5);

annotation('arrow', [.42 .64], [.35 .35],'Linewidth',1.5);
txt = annotation(fig,'textbox');
txt.String = 'Resampling'; txt.Position = [.44,.35,.17,.07];
txt.FontSize = 12; txt.EdgeColor = 'none';
txt.FaceAlpha = 0.5; txt.BackgroundColor = 'g';

savefig(fig,'Figure_supplement_6.fig','compact');