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
stepflist = [5:25];

scvs = cell(length(par),length(mnames),2);
sauc = cell(length(par),length(mnames),2);
trialnums = cell(length(par),2);
steplists = cell(length(par),2);
trialtype = cell(length(par),2);
rposall = cell(length(par),2);

trialtylist = {'correct','incorrect'};
crin = 1; dirlesit = {'bw','fw'};
for idir = 1:2
    direc = dirlesit{idir};
    for ip = [1:length(par)]
        session = par{ip,1}; trialtytmp = {'correct','incorrect'};
        resultsfile = [session,'\',trialtytmp{crin},'_SVMClassification_',direc,'.mat'];
        
        %% plotting
        load(resultsfile,'rate_cvs','rate_auc','rate_aucxy','rpos_all','TrNmids','trialty');
        for ifet = 1:4
            for ist = 1:length(stepflist)
                sauc{ip,ifet,idir} = [sauc{ip,ifet,idir},rate_auc{ist,ifet}];
                scvs{ip,ifet,idir} = [scvs{ip,ifet,idir},rate_cvs{ist,ifet}];
            end
            rposall{ip,idir} = rpos_all;
        end
        trialnums{ip,idir} = TrNmids;
        steplists{ip,idir} = stepflist;
        trialtype{ip,idir} = trialty;
    end
end

save(['Figure6_results_data_svmprediction_',trialtylist{crin}],'mnames','scvs','sauc','rposall','trialnums','trialtype','steplists');
