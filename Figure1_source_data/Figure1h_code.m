clear; clc;
warning('off');

%% load data
addpath('../Utilities');
if ~exist('Figure1h_results_data.mat','file')
    load('Figure1h_source_data.mat');

    %% Figure 1h decoding
    deme = 'vb'; recom = false;
    vidx = find(Speed>0.05);

    rtim = Timestamps(vidx);
    rpos = normalize(Position(vidx),1,'range');

    fparaw = zscore(fparw(vidx,:),[],2);
    fpadsk = zscore(fpadk(vidx,:),[],2);
    fpares = zscore(fpars(vidx,:),[],2);
    fparesN = zscore(fparsN(vidx,:),[],2);

    mnames={'raw data','despiked data','residual','normalized residual'};

    icdx = 1:128;

    %% all data
    trackLen = 290; OB=cell(0);
    OB{1} = zscore(fparaw(:,icdx));
    OB{2} = zscore(fpadsk(:,icdx));
    OB{3} = zscore(fpares(:,icdx));
    OB{4} = zscore(fparesN(:,icdx));

    % OLE
    runlist = 1:length(OB);
    for im=runlist
        obs = OB{im};
        % Generate basis...
        [basis,Xrbf] = get1Dbasis('vonmises',75,rpos*2*pi,100);

        pvec = linspace(0,2*pi,trackLen/2); % circular phase 
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
        f = zscore(f,[],2);
        idx = 1:size(f,1);
        [~,maxpost] = max(f,[],2);
        DecodedPos(:,2*(im-1)+1) = trackLen*rpos(idx);
        DecodedPos(:,2*(im-1)+2) = trackLen*maxpost/length(pvec);

        % Circular error...
        errp = circ_dist(     rpos(idx)*2*pi, maxpost/length(pvec)*2*pi)*trackLen/(2*pi);
        errm = circ_dist(2*pi-rpos(idx)*2*pi, maxpost/length(pvec)*2*pi)*trackLen/(2*pi);
        % Save results
        errpAll(:,im) = errp;
        errmAll(:,im) = errm;
        % non-directional error
        [err_nd_All(:,im),mini] = (min([abs(errmAll(:,im)) abs(errpAll(:,im))]'));
        pos_hat(:,im) = maxpost'/length(pvec)*2*pi;
    end

    % Median error in cm with s.e.
    for im=runlist
        bootstat = bootstrp(500,'median',abs(err_nd_All(:,im)));
        median_error(im,:) = [mean(bootstat),std(bootstat)];
    end
    rpos = rpos*trackLen;
    pos_hat = round(pos_hat*trackLen/(2*pi));
    save('Figure1h_results_data.mat','trackLen','rpos','pos_hat','mnames','err_nd_All','median_error');
end

%% Figure 1h
load('Figure1h_results_data.mat')

col = {[7,7,7]./255;
       [188,80,157]./255;
       [243,152,0]./255;
       [149,97,52]./255};
   
figure('Color','white');
for ih = 1:4
    h(:,1) = cdfplot(err_nd_All(:,ih)); hold on;
    set(h(:,1), 'LineStyle', '-', 'LineWidth',3, 'Color', col{ih});
    text(30,0.5-ih*0.1,sprintf('%.1f cm',median_error(ih,1)),'Color',col{ih});
end
title('figure 1h'); grid off;  xlim([0,50]);
plot([0,50],[0.5,0.5],'k:');
xlabel('Decoding error (cm)'); ylabel('CDF');
yticks([0,0.5,1]); xticks([0,25,50]);
set(gca,'LineWidth',1);

savefig(gcf,'Figure1h','compact')