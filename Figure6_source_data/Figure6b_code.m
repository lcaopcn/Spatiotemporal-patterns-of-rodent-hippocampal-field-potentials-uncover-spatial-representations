clear; clc;
warning('off');

crin = 1;
load('Figure6_results_data_decodingprediction.mat');

figure('Color','white');
for iplt = 1:2
    subplot_tight(1,2,iplt,[0.12,0.1])

    scvs_tmp = cell(4,1);
    stepflist = steplists{1,1};
    for ifet = 1:4
        for ip = [1:4]+(iplt-1)*4
            scvs_tmp{ifet} = [scvs_tmp{ifet}; squeeze(scvs{ip,ifet,crin})];
        end
    end

    h1 = plot_areaerrorbar(scvs_tmp{1},coloropt(1,'-','*',0.5)); hold on;
    h2 = plot_areaerrorbar(scvs_tmp{2},coloropt(2,'-','*',0.5)); hold on;

    yticks([0:0.1:1.0]); yticklabels([0:10:100]); ylim([0.4,1.0]);
    bkslen = length(stepflist); xlim([1,bkslen]);
    xticks([1:1:bkslen]); xtickangle(45) 
    xticklabels([['  0';num2str([stepflist(1:end-1)*0.1-0.4]','%.1f')],repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
    h = gca; h.YAxis.MinorTick = 'on';
    h.YAxis.MinorTickValues = 0.1:0.1:0.9;
    legend([h1(1),h2(1)],'FPA_{uhf}+dLFP_{\theta}','Spikes');

    xlabel('Time before choice point (s)');
    ylabel('Prediction accuracy (%)');
end

savefig(gcf,'Figure6b','compact')