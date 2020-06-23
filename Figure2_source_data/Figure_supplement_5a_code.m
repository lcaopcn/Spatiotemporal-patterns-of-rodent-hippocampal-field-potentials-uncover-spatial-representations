clear; clc;
warning('off');

example = [95,707,278]; deme = 'vb'; 

icdx = 1:128; 
OB=cell(0); rPOS=cell(0);

%% Load and PreProcess theta Data
load('Figure_supplement_5_source_data_wake.mat');
vidx = find(Speed>0.05);    % speed threshold 5 cm/s
rpos_init = normalize(Position(vidx,:),1,'range');
    
rspk = Spikes(vidx,:);
rFPAuhf = normalize(FPAuhf(vidx,:),2,'range');

OB{1,2} = rFPAuhf(:,icdx);
OB{1,1} = rspk(:,:);
rPOS{1} = rpos_init;

%% Load and PreProcess ripple Data
filterband = [140,250]; probthr = 0; recom = false;
load('Figure_supplement_5_source_data_ripple.mat');
Frvall = Frvall;
rFPAuhf = normalize(FPAuhf,2,'range');

OB{2,2} = rFPAuhf(:,icdx);
OB{2,1} = Frvall(:,:);
rPOS{2} = normalize(probpos,1,'range');

trackLen = 320;

%% OLE learn from running fpa
obs = OB{1,2};
% Generate basis...
[basis,Xrbf] = get1Dbasis('vonmises',75,rPOS{1}*2*pi,100);

w2 = zeros(0);
for j=1:size(Xrbf(:,:),2)
    [w2(:,j),V,invV,~,~,~,~,~] = ...
        VB_ARD_linear_regression(obs(:,:), Xrbf(:,j));
end

%% OLE learn from running spk
obs = OB{1,1};
% Generate basis...
[basis,Xrbf] = get1Dbasis('vonmises',75,rPOS{1}*2*pi,100);

w1 = zeros(0);
for j=1:size(Xrbf(:,:),2)
    [w1(:,j),V,invV,~,~,~,~,~] = ...
        VB_ARD_linear_regression(obs(:,:), Xrbf(:,j));
end

    
%% OLE test for ripple events
pvec = linspace(0,2*pi,trackLen/2); % circular phase 
[~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);

for ij = 1:length(example)
    ir = example(ij);
    fprintf('processing #%d (total %d) ...\n',ir,max(ReplayEventIdx));
    
    idxtr = find(ReplayEventIdx~=ir);
    idxts = find(ReplayEventIdx==ir);

    % testing
    pr{ij} = probpro(idxts,:);
    f1{ij} = OB{2,1}(idxts,:) * (w1*dbasis');
    f2{ij} = OB{2,2}(idxts,:) * (w2*dbasis');
end

%% ploting
oFig = figure('Visible','on','Color','white');
scup = 10; rows = 10; cols = 11;
pos = 0.01:0.02:3.19; spacing = [0.01,0.01]; colormap(flipud(bone));
for ij = 1:length(example)
    ir = example(ij);
    idxts = find(ReplayEventIdx==ir);
    tim = [1:length(idxts)].*0.02;

    subplot_tight(rows,cols,reshape([0:2]*cols+[2:4]',9,1)+(ij-1)*3,spacing);
    [~,puptivepos] = max(probpro(idxts,:),[],2);
    imagesc(tim,pos,smoothdata(probpro(idxts,:),2,'gaussian',10)'); hold on;
    scatter(tim,pos(puptivepos),36,'or','filled');
    set(gca,'YDir','reverse'); clim([0,0.1]);
    yticks([0,3.2]); yticklabels([]); xticks([0:0.1:0.5]); xticklabels([]);
    
    subplot_tight(rows,cols,reshape([3:5]*cols+[2:4]',9,1)+(ij-1)*3,spacing);
    [~,puptivepos] = max(f1{ij},[],2);
    sf1{ij} = exp(scup.*f1{ij})./sum(exp(scup.*f1{ij}),2);
    imagesc(tim,pos,sf1{ij}'); hold on;
    scatter(tim,pos(puptivepos),36,'sg','filled');
    set(gca,'YDir','normal'); clim([0,0.1]);
    yticks([0,3.2]); yticklabels([]);  xticks([0:0.1:0.5]); xticklabels([]);
    
    subplot_tight(rows,cols,reshape([6:8]*cols+[2:4]',9,1)+(ij-1)*3,spacing);
    [~,puptivepos] = max(f2{ij},[],2);
    sf2{ij} = exp(scup.*f2{ij})./sum(exp(scup.*f2{ij}),2);
    imagesc(tim,pos,sf2{ij}'); hold on;
    scatter(tim,pos(puptivepos),36,'^g','filled');
    set(gca,'YDir','normal'); clim([0,0.1]);
    xlabel('Time (s)'); yticks([0,3.19]); xticks([0:0.1:0.5]);
    if (ij == 1) 
        yticklabels([0,320]); yl = ylabel('Position (cm)');
        set(yl,'position', get(yl,'position')+[0.02,0,0]);
    else, yticklabels([]); end
end

% colorbar for fisrt row
cmaxrange = 0.1;
ax = axes; colormap(ax,flipud(bone)); clim([0,cmaxrange]);
cb = colorbar(ax); cb.Position = [0.06,0.53,0.023,0.3]; 
cb.Label.String = "Probability"; cb.Label.Position = [-0.6,cmaxrange/2];
set(cb,'XTick',[0,cmaxrange]); cb.FontSize = 8; ax.Visible = 'off';


savefig(oFig,'Figure_supplement_5a','compact');
