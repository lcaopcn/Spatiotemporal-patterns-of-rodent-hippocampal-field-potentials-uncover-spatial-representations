clear; clc;
warning('off');

%%
addpath('../Utilities');
if ~exist('Figure1de_results_data.mat','file')
    load('Figure1de_source_data.mat');
    xlen = 120; ylen = 120; spbin = 5; tmbin = 0.3;

    mnames={'Clustered Spike','dLFPs_{\theta}','FPA_{uhf}','dLFPs_{\theta}+FPA_{uhf}','All'};
    icdx = 1:size(FPAuhf,2);
    isk = 1:size(Spikes,2);

    vidx = Speed>mean(Speed); nfoldcv = 10;
    rpos = normalize([xPosition(vidx,:), yPosition(vidx,:)],1,'range');

    rphase = zscore(whiten(dLFPTheta_phase(vidx,:)')',[],2);
    ramp = zscore(whiten(dLFPTheta_amp(vidx,:)')',[],2);
    rFPAuhf = zscore(whiten(FPAuhf(vidx,:)')',[],2);
    rspk = zscore(Spikes(vidx,:),[],2);

    OB{1} = rspk(:,isk);
    OB{2} = [rphase(:,icdx),ramp(:,icdx)];
    OB{3} = rFPAuhf(:,icdx);
    OB{4} = [rFPAuhf(:,icdx),rphase(:,icdx),ramp(:,icdx)];
    OB{5} = [rspk(:,isk),rFPAuhf(:,icdx),rphase(:,icdx),ramp(:,icdx)];
    
    basis.n = 12;   % number of basis functions along each dimension
    [basis,Xrbf] = get2Dbasis('gaussian',[basis.n basis.n],rpos);

    % Positions to use for decoding...
    decode_gridn = xlen/spbin;
    [py,px]=meshgrid(linspace(0,1,decode_gridn),linspace(0,1,decode_gridn));
    pvec = [px(:) py(:)];
    [tmp,dbasis] = get2Dbasis('gaussian',[basis.n basis.n],pvec);

    F = cell(0); xyallhat = cell(0); errall = cell(0);

    for im = 1:length(OB)
        nfoldcv = 10;

        % Fit n-fold cross-validated, L2-regularized place-field models for each neuron/lfp-channel...
        obs = OB{im};
        for ic=1:size(obs,2)
            fprintf('\nChannel %03i/%03i... ',ic,size(obs,2))
            m(ic) = fitCVridge(Xrbf,obs(:,ic),nfoldcv,[0 logspace(-3,1,10)]);
        end
        fprintf('\n');

        % Collect place fields for each cv-fold...
        lam = cell(0);
        for i=1:nfoldcv
            for ic=1:size(obs,2)
                lam{i}(:,ic) = m(ic).breg(1,i)+dbasis*m(ic).breg(2:end,i);
            end
        end

        %% Bayesian...
        f = [];
        for i=1:nfoldcv
            % Estimate training set spike-prediction noise...
            sigma=[];
            for ic=1:size(obs,2)
                sigma(ic) = std(obs(m(1).cvidx_tr{i},ic)-(m(ic).breg(1,i)+Xrbf(m(1).cvidx_tr{i},:)*m(ic).breg(2:end,i)));
            end
            sigma(sigma<0.01)=mean(sigma);

            f(m(1).cvidx_ts{i},:) = decodeBayesian_gauss(obs(m(1).cvidx_ts{i},:),lam{i},sigma');
        end

        %% smooth kernels
        difx = pvec(1:size(f,2),1) - pvec(1:size(f,2),1)';
        dify = pvec(1:size(f,2),2) - pvec(1:size(f,2),2)';
        I = eye(2);
        
        beta = 0.06; show = false;
        kernalfile = sprintf('Guaussian_Kernel_2D_%dcm_%.2f.mat',spbin,beta);
        if ~exist(kernalfile,'file')
            if show
                figure('Position',[600,300,500,500]);
            end
            temp = zeros(0);
            for poslast = 1:size(f,2)
                diffxy = [difx(:,poslast),dify(:,poslast)]';
                temp(:,poslast) = diag(exp(-diffxy'*I*diffxy./(beta^2))'./(2*pi*beta^2));
                if show
                    % visualization
                    clf;
                    scatter(pvec(1:decode_gridn^2,1),pvec(1:decode_gridn^2,2),256,temp(:,poslast)*10,'filled');
                    xticks([]); yticks([]); title('smoothe kernel'); pbaspect([1 1 1]); drawnow;
                end
                fprintf('\nLoactioon bin: %03i/%03i... ',poslast,size(f,2));
            end
            save(kernalfile,'temp');
        else
            load(kernalfile);
        end

        f = smoothingFilter(f,temp);

        decodePostProcess_2D;
        xyallhat{im} = xyhat;
        fall{im} = f;
        errall{im} = err; 
        bootstat = bootstrp(500,'median',abs(err));
        median_error(im,1:2) = [mean(bootstat),std(bootstat)];
        err_nd_All = decoding_err;
    end
    save('Figure1de_results_data.mat','decode_gridn','median_error','err_nd_All','mnames','rpos','xyallhat');
end

%% Figure 1d
load('Figure1de_results_data.mat');
layoutall = [1,3; 2,4; 7,9; 8,10];
col = {'r','b','g','c'};
titlelist = {'Clustered spikes','dLFP_{\theta}','FPA_{uhf}','dLFP_{\theta}+FPA_{uhf}'};

figure('Color','white');
for it = 1:4
    subplot(5,2,layoutall(it,1))
    plot(rpos(:,1),'k','LineWidth',2); hold on;
    sca = scatter(1:size(xyallhat{it},1),xyallhat{it}(:,1),36,col{it}','filled');
    sca.MarkerFaceAlpha = 0.5;
    xlim([570,750]);axis off; box off;
    title(titlelist{it});

    subplot(5,2,layoutall(it,2))
    plot(rpos(:,2),'k','LineWidth',2); hold on;
    sca = scatter(1:size(xyallhat{it},1),xyallhat{it}(:,2),36,col{it},'filled');
    sca.MarkerFaceAlpha = 0.5;
    xlim([570,750]);axis off; box off;
end

savefig(gcf,'Figure1d','compact');

%% Figure 1e
col = {'r','b','g','c'};
   
figure('Color','white');
for ih = 1:4
    h(:,1) = cdfplot(err_nd_All{ih}); hold on;
    set(h(:,1), 'LineStyle', '-', 'LineWidth',3, 'Color', col{ih});
    text(30,0.5-ih*0.1,sprintf('%.1f cm',median_error(ih,1)),'Color',col{ih});
end
title('figure 1h'); grid off;  xlim([0,50]);
plot([0,50],[0.5,0.5],'k:');
xlabel('Decoding error (cm)'); ylabel('CDF');
yticks([0,0.5,1]); xticks([0,25,50]);
set(gca,'LineWidth',1);

savefig(gcf,'Figure1e','compact')