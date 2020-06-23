clear; clc;
warning('off');

%% load data
addpath('../Utilities');
if ~exist('Figure1bc_results_data.mat','file')
    load('Figure1bc_source_data.mat');

    vidx = find(Speed>0.05);    % speed threshold 5 cm/s
    rpos_init = normalize(Position(vidx,:),1,'range');

    rspk = Spikes(vidx,:);
    rphase = zscore(dLFPTheta_phase(vidx,:),[],2);
    ramp = zscore(dLFPTheta_amp(vidx,:),[],2);
    rFPAuhf = normalize(FPAuhf(vidx,:),2,'range');

    %% Figure 1b&c decoding
    deme = 'vb';
    mnames={'Spikes','dLFPs_{\theta}','FPA_{uhf}','Spikes (+history)','dLFPs_{\theta} (+history)','FPA_{uhf} (+history)'};

    median_error = zeros(0);
    DecodedPos = zeros(0);

    trackLen = 290;
    OB=cell(0); plotline = 2;
    rpos = rpos_init(2:end);
    OB{1} = zscore(rspk(2:end,:));
    OB{2} = zscore([rphase(2:end,:),ramp(2:end,:)]);
    OB{3} = zscore(rFPAuhf(2:end,:));

    OB{4} = zscore([rspk(2:end,:),rspk(1:end-1,:)]);
    OB{5} = zscore([rphase(2:end,:),ramp(2:end,:),rphase(1:end-1,:),ramp(1:end-1,:)]);
    OB{6} = zscore([rFPAuhf(2:end,:),rFPAuhf(1:end-1,:)]);

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
        errp = circ_dist(     rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
        errm = circ_dist(2*pi-rpos(idx)*2*pi, maxpost'/length(pvec)*2*pi)*trackLen/(2*pi);
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
    save('Figure1bc_results_data','trackLen','mnames','rpos','pos_hat','err_nd_All','median_error');
end

%% Figure 1b
load('Figure1bc_results_data');
plotwindow = [94.6,108];
timestamps = [1:length(rpos)]*0.1;
idx = [0;diff(rpos)] > 1.5;
rpos(idx) = nan; 
ms  = 36;

figure('Color','white');
subplot(3,1,1)
plot(timestamps,rpos,'k','LineWidth',3); hold on;
s = scatter(timestamps,pos_hat(:,1),ms,'r','filled');
s.MarkerFaceAlpha = 0.5;
xlim(plotwindow); box off; axis off;
title("Clustered spikes")

subplot(3,1,2)
plot(timestamps,rpos,'k','LineWidth',3); hold on;
s = scatter(timestamps,pos_hat(:,2),ms,'b','filled');
s.MarkerFaceAlpha = 0.5;
xlim(plotwindow); box off; axis off;
title("Demodulated theta (dLFP_{\theta})");

subplot(3,1,3)
plot(timestamps,rpos,'k','LineWidth',3); hold on;
s = scatter(timestamps,pos_hat(:,3),ms,'g','filled');
s.MarkerFaceAlpha = 0.5;
xlim(plotwindow); box off; axis off;
title(">300 Hz amplitude (FPA_{uhf})");

suptitle('Figure 1b');
savefig(gcf,'Figure1b','compact');

%% Figure 1c
figure('Color','white');

col = {[1 0 0],[0 0 1],[0 1 0]};
for it = 1:3
    isub = subplot(3,2,2*it);
    h(:,1) = cdfplot(err_nd_All(:,it)); hold on;
    h(:,2) = cdfplot(err_nd_All(:,it+3)); xlim([0,50]); title('');
    grid off; yticks([0,0.5,1]); xticks([0,25,50]);
    xlabel(''); ylabel(''); pbaspect(isub,[1 1 1]);
    if (it == 3),  xlabel('Decoding error (cm)'); ylabel('CDF'); end
    text(20,0.35,sprintf(' %.1f cm\n(%.1f cm)',median_error(it,1),median_error(it+3,1)),'Color',col{it});
    
    set(h(:,1), 'LineStyle', '-', 'LineWidth',1.5, 'Color', col{it});
    set(h(:,2), 'LineStyle', ':', 'LineWidth',1.5, 'Color', col{it});   
end

edges = linspace(0,trackLen,trackLen/2);
e{1}=edges; e{2}=edges;
ex = edges+mean(diff(edges))/2;
[mx,my] = meshgrid(ex,ex);
posbinc = 1:2:trackLen-1;

cmaxrange = 0.1;
for ib=1:3
    isub = subplot(3,2,2*ib-1);
    colormap(flipud(bone));
    [n,~]=hist3([rpos,pos_hat(:,ib)],'Edges',e);
    ns = imgaussfilt(n/2^2,2);
    imagesc(posbinc,posbinc,ns); set(gca,'YDir','normal');
    axis image; isub.CLim = [0,0.1];
    xlim([0,1]*trackLen); ylim([0,1]*trackLen);
    xticks([0,trackLen]); yticks([0,trackLen]);
    if ib == 3,  xlabel('True position (cm)'); ylabel('Estimated (cm)'); end
end
ax = axes; colormap(ax,flipud(bone)); ax.CLim = [0,0.1];
cb = colorbar(ax); cb.Position = [0.1,0.7,0.02,0.2]; 
cb.Label.String = "Density (count/cm^2)"; cb.Label.Position = [-2.5,cmaxrange/2];
set(cb,'XTick',[0,cmaxrange]); cb.FontSize = 8; ax.Visible = 'off';

suptitle('Figure 1c');
% savefig(gcf,'Figure1c','compact');
