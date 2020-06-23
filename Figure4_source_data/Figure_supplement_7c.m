clear; clc;
warning('off');

%% Figure_supplement_7c
load('Figure_supplement_7c_source_data');
vidx = find(lSpeed > 30);
rposl = lPosition(vidx,:);
lfpa = zscore(lFPAuhf(vidx,:),[],2);

vidx = find(rSpeed > 30);
rposr = rPosition(vidx,:); 
rfpa = zscore(rFPAuhf(vidx,:),[],2);

clear colmapc;  colmap = jet;
for ic = 1:3
    colmapl(:,ic) = interp1(1:size(colmap,1),colmap(:,ic),normalize(rposl,'range')*(size(colmap,1)-1)+1);
    colmapr(:,ic) = interp1(1:size(colmap,1),colmap(:,ic),normalize(rposr,'range')*(size(colmap,1)-1)+1);
end

colmapc = [colmapl;colmapr];
idxs = [zeros(size(colmapl,1),1);ones(size(colmapr,1),1)];
fpa = [lfpa; rfpa];

distanceal = 'euclidean';
nclustr = 3;

Z = tsne(fpa,'Algorithm','barneshut','Distance',distanceal,'Exaggeration',nclustr);

figure('Color','white');
subplot(3,5,[1:2,6:7,11:12])
idx = idxs == 0; 
scatter(Z(idx,1),-Z(idx,2),24,colmapc(idx,:),'filled');
set(gca, 'Color', [0 0 0]); box off; xticks([]); yticks([]);
subplot(3,5,[1:2,6:7,11:12]+3)
idx = idxs == 1;
scatter(Z(idx,1),-Z(idx,2),24,colmapc(idx,:),'filled');
set(gca, 'Color', [0 0 0]); box off; xticks([]); yticks([]);

subplot(3,5,8)
idxl = Postemp(:,4) == 1 | Postemp(:,4) == 2;
colmap = jet; posl = Postemp(idxl,1);
for ic = 1:3
    crmp(:,ic) = interp1(1:size(colmap,1),colmap(:,ic),normalize(posl,'range')*(size(colmap,1)-1)+1);
end
scatter(Postemp(idxl,3),Postemp(idxl,2),36,crmp,'filled'); hold on;

idxr = Postemp(:,4) == 1 | Postemp(:,4) == 3; post = Postemp;
post(Postemp(:,4)==3) = post(Postemp(:,4)==3) - (352-122);
colmap = jet; posr = post(idxr,1);
for ic = 1:3
    crmp(:,ic) = interp1(1:size(colmap,1),colmap(:,ic),normalize(posr,'range')*(size(colmap,1)-1)+1);
end
scatter(Postemp(idxr,3),Postemp(idxr,2),36,crmp,'filled');
box off; axis off;

savefig(gcf,'Figure_supplement_7c','compact');