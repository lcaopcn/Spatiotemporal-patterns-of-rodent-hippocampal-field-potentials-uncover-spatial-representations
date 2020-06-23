clear; clc
warning('off');

load('..\Figure4_source_data_run.mat');
rfpa = rFPAuhf;
[W,H] = nnmf(rfpa',64);
H = zscore(H,0,1);
rH = int64(resample2poission(H',10,0));
rfpa = int64(rfpa);
rspk = int64(rspk);
[~,~,bin] = histcounts(rpos,linspace(0,1,60));
for ib = 1:max(bin)
    ave(ib,:) = mean(rH(ib == bin,:));
end
[~,imax] = max(ave,[],1);
[~,ord] = sort(imax);
imagesc(ave(:,ord));
save('pyhsmm_data_run','rpos','rfpa','W','rH','rspk');

clearvars -except W;
load('..\Figure4_source_data_ripple.mat');
rfpa = rFPAuhf;
[W,H] = nnmf(rfpa',64);
H = zscore(H,0,1);
rH = int64(resample2poission(H',3,0));
rfpa = int64(rfpa);
rspk = int64(Frvall);
save('pyhsmm_data_ripple','rfpa','W','rH','rspk');
