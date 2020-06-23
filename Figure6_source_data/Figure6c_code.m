clear; clc;
warning('off');

crin = 1; trialtylist = {'correct','incorrect'};
load(['Figure6_results_data_svmprediction_',trialtylist{crin}]);

figure('Color','white');
iplist = [1:4]; idir = 1;
spacing = [0.05,0.1];

subplot_tight(2,1,1,spacing)
scvs_tmp = cell(4,1);
sauc_tmp = cell(4,1);

stepflist = steplists{1,1};
for ifet = 1:4
    for ip = iplist
        scvs_tmp{ifet} = [scvs_tmp{ifet}; squeeze(scvs{ip,ifet,crin})];
        sauc_tmp{ifet} = [sauc_tmp{ifet}; squeeze(sauc{ip,ifet,crin})];
    end
end

h1 = plot_areaerrorbar(scvs_tmp{1},coloropt(1,'-','*',0.5)); hold on;
h2 = plot_areaerrorbar(scvs_tmp{2},coloropt(2,'-','*',0.5)); hold on;

yticks([0:0.1:1.0]); yticklabels([0:10:100]); ylim([0.5,1.0]);
bkslen = length(stepflist); xlim([1,bkslen]);
xticks([1:1:bkslen]); xtickangle(45); xticklabels([]);
h = gca; h.YAxis.MinorTick = 'on';
h.YAxis.MinorTickValues = 0.1:0.1:0.9;
legend([h1(1),h2(1)],'FPA_{uhf}+dLFP_{\theta}','Spikes');
ylabel('Classification accuracy (%)');

subplot_tight(2,1,2,spacing)
rpos = [];
for ip = iplist
    rpos_tmp = [];
    for ist = 1:length(stepflist)
        rpos_tmp = [rpos_tmp,reshape(rposall{ip,idir}{ist},numel(rposall{ip,idir}{ist}),1)];
    end
    rpos = [rpos;rpos_tmp];
end
boxplot(-rpos-1.13,'orientation', 'vertical','symbol','');
set(gca,'YDir','reverse'); ylim([0,0.68]);
xticks([1:1:bkslen]); xtickangle(60); xlim([1,bkslen]+[-0.5,0.5]);
xticklabels([num2str([0,stepflist(1:end-1)*0.1-0.4]','%.1f'),...
                repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
xlabel('Time before choice point (s)');
ylabel({'Animal''s position','from choice point (m)'});
    
savefig(gcf,'Figure6c','compact')