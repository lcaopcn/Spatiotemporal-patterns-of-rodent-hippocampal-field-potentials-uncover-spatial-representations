clear; clc;
warning('off');

%% Figure_supplement_7a
load('Figure_supplement_7b_source_data');

vidx = find(Speed>min(nanmean(Speed),30));

rpos = (normalize(Position(vidx,:),1,'range'));
rmuaa = zscore(FPAuhf(vidx,:),[],2);

colmap = jet; colmapc = [];
for ic = 1:3
    colmapc(:,ic) = interp1(1:size(colmap,1),colmap(:,ic),rpos*(size(colmap,1)-1)+1);
end

idx = 1:size(rmuaa,1);
distanceal = 'euclidean';

figure('Color','white')
Z = tsne(rmuaa,'Algorithm','barneshut','Distance',distanceal,'Exaggeration',4,'NumDimensions',2);
scatter(Z(:,1),Z(:,2),24,colmapc(:,:),'filled'); hold on;
set(gca, 'Color', [0 0 0]); box off; xticks([]); yticks([]);

savefig(gcf,'Figure_supplement_7b','compact');