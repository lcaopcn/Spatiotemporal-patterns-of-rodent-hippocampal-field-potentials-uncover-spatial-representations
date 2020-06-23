clear; clc;
warning('off');

par = {'mouse1_session1';
       'mouse1_session2';
       'mouse1_session3';
       'mouse1_session4';
       
       'mouse2_session1';
       'mouse2_session2';
       'mouse2_session3';
       'mouse2_session4';
       };

nshf = 1000; 
mnames={'FPA_{uhf}+dLFP_{\theta}','Spikes','FPA_{uhf}','dLFP_{\theta}'};

scvs = cell(length(par),length(mnames),2);
ssfs = cell(length(par),length(mnames),2);
trialnums = cell(length(par),2);
steplists = cell(length(par),2);
trialtype = cell(length(par),2);

for crin = 1:2
    for ip = [1:length(par)]
        session = par{ip,1}; trialty = {'correct','incorrect'};
        resultsfile = [par{ip,1},'\',trialty{crin},'_DecodingPrediction.mat'];

        %% plotting
        load(resultsfile);
        for ifet = 1:4
            for ist = 1:length(stepflist)
                ssfs{ip,ifet,crin} = [ssfs{ip,ifet,crin},rate_shf{ist,ifet}(:,2)];
                scvs{ip,ifet,crin} = [scvs{ip,ifet,crin},rate_cvs{ist,ifet}(:,2)];
            end
        end
        trialnums{ip,crin} = TrNmids;
        steplists{ip,crin} = stepflist;
        trialtype{ip,crin} = trialty;
    end
end

save('Figure6_results_data_decodingprediction','mnames','scvs','ssfs','trialnums','trialtype','steplists');
