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
for crin = 1:2 % 1, correct trials; 2, incorrect trials
stepflist = [5:15];

for ip = [5]
    rate = cell(length(stepflist),4);
    rate_cvs = cell(length(stepflist),4);
    rate_shf = cell(length(stepflist),4);
    
    session = par{ip,1}; trialtype = {'correct','incorrect'};
    resultsfile = [par{ip,1},'\',trialtype{crin},'_DecodingPrediction.mat'];
    if ~exist(resultsfile,'file')
        for istepf = 1:length(stepflist)
            clearvars -except crin par ip istepf stepflist nshf ...
                    rate rate_cvs rate_shf resultsfile session trialtype;
            intansr = 20000; nchan = 64; warning('off');
            disp(['running session: ',session,' ...']);
            tic; warning('off');

            marker = [1,0]; 

            stepf = stepflist(istepf);
            weight = [1,0.9,0.8,0.7,0.6,0.5,0.5,0.5,0.4,0.3,0.2,0.1];
            dt = 0.1; deme = 'ls'; 

            %% Load and PreProcess Data
            load([session,'\Features_data'],'sidearm','centerarm','dcs','TrNm',...
                                'Spike','dLFPTheta_amp','dLFPTheta_phase','FPAuhf');
            load(sprintf('%s\\Label_data_stepf_%d',session,stepf),'pos');


            trackLen = 2.5;
            grdth = centerarm(centerarm(:,6)==marker(crin),5);
            rdcs = dcs{crin}(:,:);
            rspk = Spike{crin}(:,:);
            rpos = (pos{crin}-min(pos{crin}))/trackLen;
            rphase = zscore(dLFPTheta_amp{crin}(:,:),[],2);
            ramp = zscore(dLFPTheta_phase{crin}(:,:),[],2);
            rfpa = zscore(FPAuhf{crin}(:,:),[],2);

            %% decoding  all data OLE
            mnames={'FPA_{uhf}+dLFP_{\theta}','Spikes','FPA_{uhf}','FPA_{uhf}+dLFP_{\theta}'};
            icdx = 1:nchan;
            OB=cell(0);
            OB{1} = zscore([rfpa(:,icdx),rphase(:,icdx),ramp(:,icdx)]);
            OB{2} = zscore([rspk(:,:)]);
            OB{3} = zscore([rfpa(:,icdx)]);
            OB{4} = zscore([rphase(:,icdx),ramp(:,icdx)]);

            runlist = 1:length(OB);
            for im=runlist
                obs = OB{im};
                % Generate basis...
                [basis,Xrbf] = get1Dbasis('vonmises',75,rpos*2*pi,400);

                pvec = linspace(0,2*pi,trackLen/0.02); % circular phase 
                [~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);

                % Cross-validated predictions
                mcvp.cv = getCVidx(size(obs,1),10); f=[];
                for it=1:mcvp.cv.nfoldcv
                    if strcmp(deme,'ls')
                        % least squared estimation  training
                        w = obs(mcvp.cv.tr{it},:) \ Xrbf(mcvp.cv.tr{it},:);
                    elseif strcmp(deme,'vb')
                        % Variational Bayes training
                        w = zeros(0);
                        for j=1:size(Xrbf(mcvp.cv.tr{it},:),2)
                            [w(:,j),V,~,~,~,~,~,~] = ...
                                VB_ARD_linear_regression(obs(mcvp.cv.tr{it},:), Xrbf(mcvp.cv.tr{it},j));
                        end
                    end
                    % testing
                    f(mcvp.cv.ts{it},:) = obs(mcvp.cv.ts{it},:) * (w*dbasis');
                end

                % Plotting and error calcs... 
                f = zscore(f,[],2); idx=1:size(f,1);
                [~,maxpost]=max(f,[],2);

                [rate{istepf,im}(1,1:2),rate{istepf,im}(1,3:4)] = accuracy(maxpost,TrNm{crin},stepf,grdth,weight);
                printf('No %d [%d] decoding finished, starting cvs: ',istepf,im);
                
                
                cvs = zeros(nshf,4); wt = cell(nshf,1); mcvpcv = cell(nshf,1);
                % Shuffle cvs
                gcp(); pctRunOnAll warning('off');
                parfor ish = 1:nshf
                    % Cross-validated predictions
                    mcvpcv{ish} = getCVidx(size(obs,1),10,true); f=[];
                    for it=1:mcvpcv{ish}.nfoldcv
                        if strcmp(deme,'ls')
                            % least squared estimation  training
                            wt{ish} = obs(mcvpcv{ish}.tr{it},:) \ Xrbf(mcvpcv{ish}.tr{it},:);
                        elseif strcmp(deme,'vb')
                            % Variational Bayes training
                            wt{ish} = zeros(0);
                            for j=1:size(Xrbf(mcvpcv{ish}.tr{it},:),2)
                                [wt{ish}(:,j),V,~,~,~,~,~,~] = ...
                                    VB_ARD_linear_regression(obs(mcvpcv{ish}.tr{it},:), Xrbf(mcvpcv{ish}.tr{it},j));
                            end
                        end
                        % testing
                        f(mcvpcv{ish}.ts{it},:) = obs(mcvpcv{ish}.ts{it},:) * (wt{ish}*dbasis');
                    end

                    % Plotting and error calcs... 
                    f = zscore(f,[],2); idx=1:size(f,1);
                    [~,maxpost]=max(f,[],2);

                    cvstmp = zeros(1,4);
                    [cvstmp(1:2),cvstmp(3:4)] = accuracy(maxpost,TrNm{crin},stepf,grdth,weight);
                    cvs(ish,:) = cvstmp;
                    %printf('No %d [%d] cvs %d finished (totally %d).',istepf,im,ish,nshf);
                end
                printf('No %d [%d] cvs %d finished (totally %d).',istepf,im,nshf,nshf);
                rate_cvs{istepf,im} = cvs;

                % Shuffle
                shf = zeros(nshf,4);
                for ish = 1:nshf
                    % Cross-validated predictions
                    f=[];
                    for it=1:mcvpcv{ish}.nfoldcv
                        wdb = wt{ish}*dbasis';
                        odf = randperm(size(wdb,1));
                        wdb = wdb(odf,:);
                        % testing
                        f(mcvpcv{ish}.ts{it},:) = obs(mcvpcv{ish}.ts{it},:) * (wdb);
                    end

                    % Plotting and error calcs... 
                    f = zscore(f,[],2); idx=1:size(f,1);
                    [~,maxpost]=max(f,[],2);

                    shftmp = zeros(1,4);
                    [shftmp(1:2),shftmp(3:4)] = accuracy(maxpost,TrNm{crin},stepf,grdth,weight);
                    shf(ish,:) = shftmp;
                    %printf('No %d [%d] shffle %d finished (totally %d).',istepf,im,ish,nshf);
                end
                printf('No %d [%d] shffle %d finished (totally %d).',istepf,im,nshf,nshf);
                rate_shf{istepf,im} = shf;
            end
            toc;
        end
        TrNmids = TrNm{crin};
        trialty = trialtype{crin};
        save(resultsfile,'rate','rate_cvs','rate_shf','stepflist','TrNmids','trialty');
    end
    
    %% plotting
    load(resultsfile,'rate','rate_cvs','rate_shf','stepflist','TrNmids','trialty');
    modlist = {'u','d'};
     
    for imod = 1
        mod = modlist{imod};
        fpadcscv = []; fpadcssf = [];
        lfpdcscv = []; lfpdcssf = []; 
        spkdcscv = []; spkdcssf = [];
        for ist = 1:length(stepflist)
            fpadcssf = [fpadcssf,rate_shf{ist,3}(:,imod)];
            lfpdcssf = [lfpdcssf,rate_shf{ist,4}(:,imod)];
            spkdcssf = [spkdcssf,rate_shf{ist,2}(:,imod)];
            fpadcscv = [fpadcscv,rate_cvs{ist,3}(:,imod)];
            lfpdcscv = [lfpdcscv,rate_cvs{ist,4}(:,imod)];
            spkdcscv = [spkdcscv,rate_cvs{ist,2}(:,imod)];
        end
        
        figure('Color','white');

        h4 = plot_areaerrorbar(fpadcssf,coloropt(3,'--','none',0.3)); hold on;
        h5 = plot_areaerrorbar(lfpdcssf,coloropt(1,'--','none',0.3)); hold on;
        h6 = plot_areaerrorbar(spkdcssf,coloropt(2,'--','none',0.3)); hold on;
        
        h1 = plot_areaerrorbar(fpadcscv,coloropt(3,'-','none',0.5)); hold on;
        h2 = plot_areaerrorbar(lfpdcscv,coloropt(1,'-','none',0.5)); hold on;
        h3 = plot_areaerrorbar(spkdcscv,coloropt(2,'-','none',0.5)); hold on;
        
        yticks([0.1:0.2:0.9]); yticklabels([10:20:90]); box on;
        bkslen = length(stepflist);
        xticks([1:1:bkslen]); xtickangle(45) 
        xticklabels([num2str([0,stepflist(1:end-1)*0.1-0.4]','%.1f'),repmat('-',bkslen,1),num2str(stepflist(1:end)'*0.1,'%.1f')]);
        h = gca; h.YAxis.MinorTick = 'on';
        h.YAxis.MinorTickValues = 0.1:0.1:0.9;
        legend([h1(1),h2(1),h3(1),h4(1),h5(1),h6(1)],'FPA_{uhf}','dLFP_{\theta}','Spikes',...
                                'FPA_{uhf} shuffle','dLFP_{\theta} Shuffle','Spikes Shuffle');

        xlabel('Time before choice point (s)');
        ylabel('Prediction accuracy (%)');
        title(sprintf('%d %s trials',length(unique(TrNmids)),trialty));

        savefig(gcf,[session,'\',trialtype{crin},'_DecodingPrediction_Joint_',mod,'.fig'],'compact');
    end
    close all;
end

end

%% function
function [acc,prdrate] = accuracy(maxpost,TrNm,stepf,grdth,weight)
    %% up (1, last 5 predicted) down (2, first 5 predicted)
    len = 5;
    stalist = [len, stepf];
    for imod = 1:2
        sta = stalist(imod);
        predict = zeros(0);
        Sleft = zeros(0);
        Sright = zeros(0);
        trialsnum = unique(TrNm);
        for it = 1:length(trialsnum)
            ids = find(TrNm == trialsnum(it));
            predictpos = maxpost(ids(end-sta+[1:len]));
            Sleft(it,1) = sum((predictpos > 110) .* weight(len:-1:1)');
            Sright(it,1) = sum((predictpos > 60 & predictpos < 80) .* weight(len:-1:1)');

            if (Sleft(it,1) > Sright(it,1))
                predict(it,1) = 0;
            elseif (Sleft(it,1) < Sright(it,1)) 
                predict(it,1) = 1;
            else
                predict(it,1) = -2;
            end
        end
        acc(1,imod) = sum(predict>=0 & predict==grdth(:,1))/length(grdth);
        prdrate(1,imod) = sum(predict<0)/length(grdth);
    end
end