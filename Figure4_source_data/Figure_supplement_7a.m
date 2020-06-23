clear; clc;
warning('off');

%% Figure_supplement_7a
load('Figure_supplement_7a_source_data');

vidx = find(Speed>min(nanmean(Speed),30)); 

rpos = (normalize(Position(vidx,:),1,'range'));
rmua = zscore(FPAuhf(vidx,:),[],2);

lidx = rpos>0.5;
ridx = rpos<0.5;

colmap = jet; colmapc = [];
for ic = 1:3
    colmapc(:,ic) = interp1(1:size(colmap,1),colmap(:,ic),rpos*(size(colmap,1)-1)+1);
end

idx = 1:size(rmua,1);
distanceal = 'euclidean';

Z = tsne(rmua,'Algorithm','barneshut','Distance',distanceal,'Exaggeration',4,'NumDimensions',2);
scatter(Z(ridx,1),Z(ridx,2),24,colmapc(ridx,:),'filled'); hold on;
scatter(Z(lidx,1),Z(lidx,2),24,colmapc(lidx,:),'filled');
set(gca, 'Color', [0 0 0]); box off; xticks([]); yticks([]);

savefig(gcf,'Figure_supplement_7a','compact');