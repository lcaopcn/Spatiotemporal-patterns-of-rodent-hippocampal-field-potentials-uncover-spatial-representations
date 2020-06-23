clear; clc;
warning('off');

%% load data
addpath('../Utilities');
if ~exist('Figure5d_results_data.mat','file')
    load('Figure5d_source_data_a.mat');
    load('Figure5d_source_data_b.mat');
    deme = 'vb';
    vidx = find(Speed>5 & Speed<200); % speed threshold 30cm/s
    trialNm = trialNm(vidx);
    
    rpos = normalize(Position(vidx,:),1,'range');
    rspk = Spike(vidx,:);

    rTheta_dphase = zscore(Theta_dphase(vidx,:),[],2);
    rTheta_damp = zscore(Theta_damp(vidx,:),[],2);
    
    rTheta_dmdu = [rTheta_dphase,rTheta_damp];
    rTheta_ofpat = zscore(Theta_ofpat(vidx,:),[],2);
    
    rFPAuhf = zscore(FPAuhf(vidx,:),[],2);
    rUhf_ofpat = zscore(Uhf_ofpat(vidx,:),[],2);

    %% OLE decoding data
    trackLen = 580;
    OB=cell(0); plotline = 2;
    rpos = normalize(rpos,'range');
    OB{1} = zscore([rTheta_dmdu]);
    OB{2} = zscore([rTheta_ofpat]);
    OB{3} = zscore([rFPAuhf]);
    OB{4} = zscore([rUhf_ofpat]);

    OB{5} = zscore([rTheta_dmdu,rFPAuhf]);
    OB{6} = zscore([rTheta_ofpat,rUhf_ofpat]);
    OB{7} = zscore([rTheta_dmdu,rTheta_ofpat]);
    OB{8} = zscore([rFPAuhf,rUhf_ofpat]);

    mnames = {'1','2','3','4','1+3','2+4','1+2','3+4','1+2+3+4'};

    runlist = 1:length(OB);
    for im=runlist
        obs = OB{im};
        % Generate basis...
        [basis,Xrbf] = get1Dbasis('vonmises',75,rpos*2*pi,100);

        pvec = linspace(0,2*pi,trackLen/5); % circular phase 
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
        f = zscore(f')';
        idx=1:size(f,1);
        [~,maxpost]=max(f');
        DecodedPos(:,2*(im-1)+1) = trackLen*rpos(idx);
        DecodedPos(:,2*(im-1)+2) = trackLen*maxpost/length(pvec);

        % Circular error...
        err =  circ_dist(     rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
        errp = circ_dist(     rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
        errm = circ_dist(2*pi-rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
        % Save results
        errAll(:,im) = err;
        errpAll(:,im) = errp;
        errmAll(:,im) = errm;
        % non-directional error
        [err_nd_All(:,im),mini] = (min([abs(errmAll(:,im)) abs(errpAll(:,im))]'));
        pos_hat(:,im) = maxpost'/length(pvec)*2*pi;
    end
    close all;

    % Median error in cm with s.e. & density heat map
    for im=runlist
        bootstat = bootstrp(500,'median',abs(err_nd_All(:,im)));
        median_error(im,:) = [mean(bootstat),std(bootstat)];
    end
    rpos = rpos*trackLen;
    pos_hat = round(pos_hat*trackLen/(2*pi));
    save('Figure5d_results_data.mat','rpos','pos_hat','mnames','DecodedPos','median_error','err_nd_All','trackLen','trialNm');
end

%% Figure 5d
load('Figure5d_results_data.mat');

figure('Color','white');
b = bar([1:8],median_error(1:8,1),'EdgeColor','none'); hold on;
xticks([1:8]); xticklabels(mnames);
erb = errorbar([1:8],median_error(1:8,1),median_error(1:8,2));
erb.LineStyle = 'none'; erb.LineWidth = 2;
erb.Color = [0,0,0.5]; ylim([4,6]); yticks([4:0.5:6]);
xlim([0.5,8.5]); box on;
xtickangle(0); ylim([3.5,5.5])
ylabel({'Median decoding error (cm)'})

savefig(gcf,'Figure5d','compact')