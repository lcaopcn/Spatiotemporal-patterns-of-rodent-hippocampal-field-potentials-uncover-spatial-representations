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
dirlesit = {'bw','fw'};

for idir = 1:2
    direc = dirlesit{idir};
    
for crin = 1 % 1, correct trials; 2, incorrect trials
   
len = 5;
for ip = [1:8]
    %% figure
    stepflist = [len:25];
    rate = cell(length(stepflist),4);
    rate_cvs = cell(length(stepflist),4);
    rate_auc = cell(length(stepflist),4);
    rate_aucxy = cell(length(stepflist),4);
    rpos_all = cell(length(stepflist),1);

    session = par{ip,1}; trialtype = {'correct','incorrect'};
    resultsfile = [session,'\',trialtype{crin},'_SVMClassification_',direc,'.mat'];
    if ~exist(resultsfile,'file')
        for istepf = 1:length(stepflist)
            clearvars -except crin direc len par ip istepf stepflist nshf ...
                rate rate_cvs rate_auc rate_aucxy rpos_all resultsfile session trialtype;
            disp(['running session: ',session,' ...']);
            tic;
            trialtype = {'correct','incorrect'};
            marker = [1,0]; 

            stepf = stepflist(istepf);

            %% Load and Process Data
            if strcmp(direc,'fw')
                load([session,'\Features_data_fw'],'sidearm','centerarm','post','dcs','TrNm',...
                            'Spike','dLFPTheta_amp','dLFPTheta_phase','FPAuhf');
            else
                load([session,'\Features_data'],'sidearm','centerarm','post','dcs','TrNm',...
                                'Spike','dLFPTheta_amp','dLFPTheta_phase','FPAuhf');
            end

            trackLen = 2.5; rpos = [];
            rdcs = []; rspk = []; rmua = []; ramp = []; rphase = [];
            grdth = centerarm(centerarm(:,6)==marker(crin),5);
            trialsnum = unique(TrNm{crin});
            for it = 1:length(trialsnum)
                ids = find(TrNm{crin} == trialsnum(it));
                if length(ids) < stepf, continue; end
                switch direc
                    case 'fw'
                        idxs = ids([1:len]+stepf-len);
                    case 'bw'
                        idxs = ids(end-stepf+[1:len]);
                end
                rdcs = [rdcs;unique(dcs{crin}(idxs,:))];
                rpos = [rpos;reshape(post{crin}(idxs,:),1,numel(post{crin}(idxs,:)))];
                rspk = [rspk;reshape(zscore(Spike{crin}(idxs,:),[],2),1,numel(Spike{crin}(idxs,:)))];
                rmua = [rmua;reshape(zscore(FPAuhf{crin}(idxs,:),[],2),1,numel(FPAuhf{crin}(idxs,:)))];
                ramp = [ramp;reshape(zscore(dLFPTheta_amp{crin}(idxs,:),[],2),1,numel(dLFPTheta_amp{crin}(idxs,:)))];
                rphase = [rphase;reshape(zscore(dLFPTheta_phase{crin}(idxs,:),[],2),1,numel(dLFPTheta_phase{crin}(idxs,:)))];
            end
            rpos_all{istepf,1} = rpos;

            mnames={'FPA_{uhf}+dLFP_{\theta}','Spikes','FPA_{uhf}','dLFP_{\theta}'};
            OB=cell(0);
            OB{1} = zscore([rmua(:,:),rphase(:,:),ramp(:,:)]);
            OB{2} = zscore([rspk(:,:)]);
            OB{3} = zscore([rmua(:,:)]);
            OB{4} = zscore([rphase(:,:),ramp(:,:)]);

            runlist = 1:length(OB);
            for im=runlist
                obs = OB{im};

                cvs = zeros(nshf,1);
                % Shuffle cvs
                prdit = []; auroc = []; xauc = []; yauc = [];
                parfor ish = 1:nshf
                    score = [];
                    % Cross-validated predictions
                    mcvpcv = getCVidx(size(obs,1),10,true); prdit=[];
                    for it=1:mcvpcv.nfoldcv
                        % training
                        SVMModel = fitcsvm(obs(mcvpcv.tr{it},:),rdcs(mcvpcv.tr{it}),...
                                    'KernelFunction','linear','Standardize',true);
                        % testing
                        [prdit(mcvpcv.ts{it},:),score(mcvpcv.ts{it},:)] = predict(SVMModel,obs(mcvpcv.ts{it},:));
                    end
                    [xauc,yauc,~,auroc(ish,1)] = perfcurve(rdcs,score(:,2),1);
                    cvs(ish) = sum(prdit == rdcs)/length(prdit);
                end
                rate_cvs{istepf,im} = cvs;
                rate_auc{istepf,im} = auroc;
                rate_aucxy{istepf,im} = [xauc,yauc];
                printf('No %d [%d] cvs %d finished (totally %d).',istepf,im,nshf,nshf);
            end
            toc;
        end
        TrNmids = TrNm{crin};
        trialty = trialtype{crin};
        save(resultsfile,'rate_cvs','rate_auc','rate_aucxy','rpos_all','TrNmids','trialty');
    end
    %% plotting
    load(resultsfile);

    fpadcscv = []; lfpdcscv = [];  spkdcscv = [];
    rpos = [];
    for ist = 1:length(stepflist)
        fpadcscv = [fpadcscv,rate_cvs{ist,1}(:,1)];
        lfpdcscv = [lfpdcscv,rate_cvs{ist,2}(:,1)];
        spkdcscv = [spkdcscv,rate_cvs{ist,3}(:,1)];
        if ist>2
            rpos_tmp = reshape(rpos_all{ist},1,numel(rpos_all{ist}));
            if numel(rpos_all{ist}) < numel(rpos_all{1})
                rpos_tmp(numel(rpos_all{ist})+1:numel(rpos_all{1})) = -1.2;
            end
        else
            rpos_tmp = reshape(rpos_all{ist},1,numel(rpos_all{ist}));
        end
        rpos = [rpos;rpos_tmp];
    end

    figure('Color','white');
    subplot(2,1,1)
    h1 = plot_areaerrorbar(fpadcscv,coloropt(3,'-','none',0.5)); hold on;
    h2 = plot_areaerrorbar(lfpdcscv,coloropt(1,'-','none',0.5)); hold on;
    h3 = plot_areaerrorbar(spkdcscv,coloropt(2,'-','none',0.5)); hold on;

    yticks([0.1:0.1:1]); yticklabels([10:10:100]); box on; ylim([0.3,1]); 
    bkslen = length(stepflist);
    xticks([1:1:bkslen]); xtickangle(60); xlim([1,bkslen]+[-0.5,0.5]);
    xticklabels([num2str([0,stepflist(1:end-1)*0.1-0.4]','%.1f'),...
                    repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
    h = gca; % Get axis to modify
    h.YAxis.MinorTick = 'on'; % Must turn on minor ticks if they are off
    h.YAxis.MinorTickValues = 0.1:0.1:0.9; % Minor ticks which don't line up with majors
    legend([h1(1),h2(1),h3(1)],'FPA_{uhf}','dLFP_{\theta}','Spikes','Location','best');
    switch direc
        case 'fw'
            xlabel('Time after starting point (s)');
        case 'bw'
            xlabel('Time before choice point (s)');
    end
    ylabel('Classification accuracy (%)');
    title(sprintf('%d %s trials',length(unique(TrNmids)),trialty));

    subplot(2,1,2)
    boxplot(-rpos'-1.15,'orientation', 'vertical','symbol','');
    set(gca,'YDir','reverse'); ylim([0,0.65]);
    xticks([1:1:bkslen]); xtickangle(60); xlim([1,bkslen]+[-0.5,0.5]);
    xticklabels([num2str([0,stepflist(1:end-1)*0.1-0.4]','%.1f'),...
                    repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
    switch direc
        case 'fw'
            xlabel('Time after starting point (s)');
        case 'bw'
            xlabel('Time before choice point (s)');
    end
    ylabel({'Animal''s position','from choice point (m)'});

    savefig(gcf,[session,'\',trialtype{crin},'_SVMClassification_pos_',direc,'.fig'],'compact');
    close all;
end

end

end