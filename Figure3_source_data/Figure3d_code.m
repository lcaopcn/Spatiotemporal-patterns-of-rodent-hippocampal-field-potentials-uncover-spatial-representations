clear; clc;
warning("off");

%% load data
raw_ripple = load("Figure3d_source_data_raw_ripple.txt");
timestamps = [1:length(raw_ripple)]'/1250;

decodedpos = load("Figure3d_source_data_out_OLE.txt");
binstamps =  [1:size(decodedpos,1)]'/50;
epochs = load("Figure3d_source_data_raw_epochs.txt");

eval = load("Figure3d_source_data_out_Evaluation.txt");
rippleAmp = eval(:,2);
evalue = eval(:,3);
pval = eval(:,4);
score = eval(:,5);
detect = eval(:,6);
runtime = eval(:,7);
decodedpos = exp(decodedpos)./sum(exp(decodedpos),2);

ie = 9;
idx = [1:3125]+(ie-1)*3125;
timerange = timestamps(idx);
findtime = binstamps(find(score>13.8115 & binstamps>timerange(1) & binstamps<timerange(end),1,'first'));
idxbin = [1:125]+(ie-1)*125;
binsrange = binstamps(idxbin)-binstamps(idxbin(2));

evalue_tmp = evalue(idxbin);
detect_tmp = detect(idxbin);
idxsta = find(diff(evalue_tmp)==1);
idxsto = find(diff(evalue_tmp)==-1);
if ~isempty(idxsta) && idxsta(1) < 3
    idxsta(1) = [];
    idxsto(1) = [];
end
sta = binsrange(idxsta-2);
sto = binsrange(idxsto+1);
numevt = min(length(sta),length(sto));

idxrang = [1:floor(epochs(ie)/0.02)]+50;
candidate(ie) = length(find(evalue_tmp(idxrang)));
significant(ie) = length(find(detect_tmp(idxrang)));


ripthr = 99.7942;
%% plot
fig = figure("Color",'w');
fig.OuterPosition = [160,80,1620,920];
fig.InnerPosition = [160,80,1600,920];

subplot(5,1,1)
plot(timerange,zscore(double(raw_ripple(idx))),'k');
ylim([-4,4]); yticks([-3,0,3]); xlim(timerange(1)+[0.02,2]);
ylabel({'Raw voltage','(z-scored)'}); xticks([]);

subplot(5,1,2)
plot(binsrange,rippleAmp(idxbin),'b'); hold on;
plot(binsrange([1,length(binsrange)]),[1,1]*ripthr,'b:');
xlim(binsrange([1,length(binsrange)])); xlim([-0.02,1.96]);
ylabel({'Ripple Band','Amplitude'}); xticks([]);

subplot(5,1,3)
imagesc(binstamps(idxbin),[1:size(decodedpos(idxbin,:),2)]*0.05,decodedpos(idxbin,:)');
hold on; colormap(flipud(bone));
if ~isempty(findtime)
    plot([findtime,findtime],[0,2.9],'LineWidth',2,'Color','g');
end
set(gca,'YDir','normal'); xlim(timerange(1)+[0.02,2]);
ylabel('Position (cm)'); xticks([]); yticks([0.05,1.5,2.9]); yticklabels([0,150,290])
                
subplot(5,1,4);
yyaxis right
pval(pval==1) = nan; pval(pval<1e-4) = 1e-5;
plot(binstamps(idxbin([1,length(idxbin)])),[0.05,0.05],'r:','LineWidth',2); hold on;
plot(binstamps(idxbin),pval(idxbin),'r-','LineWidth',2); ylim([1e-3,2]);
set(gca,'Yscale','log'); yticks([1e-3,1e-2,1e-1,1e0])
ylabel('P-value');
yyaxis left;
plot(binstamps(idxbin),score(idxbin),'b','LineWidth',2); hold on;
plot(binstamps(idxbin([1,length(idxbin)])),[-3*log(0.01),-3*log(0.01)],'b:','LineWidth',2); hold on;
if ~isempty(findtime)
    plot([findtime,findtime],[0,15],'LineWidth',2,'Color','g');
end
xlim(timerange(1)+[0.02,2]); ylabel('Score'); xticks([]);
ylim([0,15]);

subplot(5,1,5);
bar(binstamps(idxbin),runtime(idxbin),'EdgeColor','none'); hold on;
plot(binstamps(idxbin([1,length(idxbin)])),[1,1]*20,'k:','LineWidth',2);
ylim([0,25]);
xlim(timerange(1)+[0.02,2]);
xticks(timerange(1)+0.1+[0,0.5,1.0,1.5,2.0]); xticklabels([0,0.5,1.0,1.5,2.0]);
xlabel('Time (s)'); ylabel({'Computation','Time (ms)'}); 

savefig(gcf,'Figure3d','compact');
