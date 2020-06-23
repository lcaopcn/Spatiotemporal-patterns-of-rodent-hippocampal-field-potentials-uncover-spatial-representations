clear; clc;
warning('off');

crin = 1; trialtylist = {'correct','incorrect'};
load(['Figure6_results_data_svmprediction_',trialtylist{crin}]);
anl = {'a','b'};
stepflist = [5:25];

for ia = 1:2
figure('Color','white','Position',[0,0,1600,600]);
direstr = {'Backward','Forward'};
dirlesit = {'bw','fw'};
for idir = 1:2
    direc = dirlesit{idir};
    for ip = [1:4]+(ia-1)*4
        subplot_tight(2,4,ip-(ia-1)*4+(idir-1)*4,[0.15,0.05])

        scvs_tmp = cell(4,1);
        for ifet = 1:4
            scvs_tmp{ifet} = squeeze(scvs{ip,ifet,idir});
            peak(ip,ifet,idir) = mean(scvs_tmp{ifet}(:,1));
        end
        
        h1 = plot_areaerrorbar(scvs_tmp{3},coloropt(3,'-','*',0.5)); hold on;
        h2 = plot_areaerrorbar(scvs_tmp{4},coloropt(1,'-','*',0.5)); hold on;
        h3 = plot_areaerrorbar(scvs_tmp{2},coloropt(2,'-','*',0.5)); hold on;

        yticks([0:0.2:1.0]); yticklabels([0:20:100]); ylim([0.3,1.0]);
        bkslen = length(stepflist); xlim([1,bkslen]);
        xticks([1:1:bkslen]); xtickangle(90); xticklabels([]);
        if idir == 1||2
            xticklabels([['  0';num2str([stepflist(1:end-1)*0.1-0.4]','%.1f')],...
                        repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
        end
        if idir == 1
            xlabel('Time before choice point (s)');
        elseif idir == 2
            xlabel('Time after start point (s)');
        end
        h = gca; h.YAxis.MinorTick = 'on';
        h.YAxis.MinorTickValues = 0.1:0.1:0.9;
        if ip == 1 && idir == 1
            legend([h1(1),h2(1),h3(1)],'FPA_{uhf}','dLFP_{\theta}','Spikes');
        end
        
        ylabel('Prediction accuracy (%)');
        title(sprintf('session %d, %d %s trials, %s',ip-(ia-1)*4,length(unique(trialnums{ip,crin})),trialtype{ip,idir},direstr{idir}));
    end
end
savefig(gcf,['Figure_supplement_10',anl{ia}],'compact');
end