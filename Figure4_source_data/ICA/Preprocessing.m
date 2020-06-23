clear; clc
warning('off');

load('..\Figure4_source_data_run.mat');
rfpa = rFPAuhf;
ri = rica(rfpa,64);
H = rfpa*ri.TransformWeights;
H = zscore(H,0,1);
rH = int64(resample2poission(H,5,0));
rfpa = int64(rfpa);
rspk = int64(rspk);
[~,~,bin] = histcounts(rpos,linspace(0,1,60));
for ib = 1:max(bin)
    ave(ib,:) = mean(rH(ib == bin,:));
end
[~,imax] = max(ave,[],1);
[~,ord] = sort(imax);
imagesc(ave(:,ord));
save('pyhsmm_data_run','rpos','rfpa','ri','rH','rspk');

clear;
load('..\Figure4_source_data_ripple.mat');
rfpa = rFPAuhf;
ri = rica(rfpa,64);
H = rfpa*ri.TransformWeights;
H = zscore(H,0,1);
rH = int64(resample2poission(H,5,0));
rfpa = int64(rfpa);
rspk = int64(Frvall);
save('pyhsmm_data_ripple','rfpa','ri','rH','rspk');
