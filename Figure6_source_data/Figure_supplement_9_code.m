clear; clc;
warning('off');

load('Figure6_results_data_decodingprediction.mat');
anl = {'a','b'};

for ia = 1:2
figure('Color','white','Position',[0,0,1600,600]);
for crin = 1:2
    for ip = [1:4]+(ia-1)*4
        subplot_tight(2,4,ip-(ia-1)*4+(crin-1)*4,[0.15,0.05])

        scvs_tmp = cell(4,1);
        ssfs_tmp = cell(4,1);
        stepflist = steplists{1,1};
        for ifet = 1:4
            scvs_tmp{ifet} = squeeze(scvs{ip,ifet,crin});
            ssfs_tmp{ifet} = squeeze(ssfs{ip,ifet,crin});
            peak(ip,ifet,crin) = mean(scvs_tmp{ifet}(:,1));
        end
        
        h4 = plot_areaerrorbar(ssfs_tmp{3},coloropt(3,'--','none',0.5)); hold on;
        h5 = plot_areaerrorbar(ssfs_tmp{4},coloropt(1,'--','none',0.5)); hold on;
        h6 = plot_areaerrorbar(ssfs_tmp{2},coloropt(2,'--','none',0.5)); hold on;
        
        h1 = plot_areaerrorbar(scvs_tmp{3},coloropt(3,'-','*',0.5)); hold on;
        h2 = plot_areaerrorbar(scvs_tmp{4},coloropt(1,'-','*',0.5)); hold on;
        h3 = plot_areaerrorbar(scvs_tmp{2},coloropt(2,'-','*',0.5)); hold on;

        yticks([0:0.2:1.0]); yticklabels([0:20:100]); ylim([0,1.0]);
        bkslen = length(stepflist); xlim([1,bkslen]);
        xticks([1:1:bkslen]); xtickangle(45) 
        xticklabels([['  0';num2str([stepflist(1:end-1)*0.1-0.4]','%.1f')],repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
        h = gca; h.YAxis.MinorTick = 'on';
        h.YAxis.MinorTickValues = 0.1:0.1:0.9;
        if ip == 1 && crin == 1
            legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],...
                        'FPA_{uhf}','dLFP_{\theta}','Spikes',...
                        'FPA_{uhf} shuffle','dLFP_{\theta} shuffle','Spikes shuffle');
        end
        xlabel('Time before choice point (s)');
        ylabel('Prediction accuracy (%)');
        title(sprintf('session %d, %d %s trials',ip-(ia-1)*4,length(unique(trialnums{ip,crin})),trialtype{ip,crin}));
    end
end
savefig(gcf,['Figure_supplement_9',anl{ia}],'compact');
end