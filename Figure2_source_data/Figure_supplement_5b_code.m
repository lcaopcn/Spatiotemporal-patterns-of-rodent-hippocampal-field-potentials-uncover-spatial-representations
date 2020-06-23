clear; clc;
warning("off");

%% load data
load('Figure_supplement_5_results_data_replayscore.mat','wZs');

wsZs(:,1) = max(wZs{1}(:,1),[],2);  
wsZs(:,2) = max(wZs{2},[],2);  
wsZs(:,3) = max(wZs{3},[],2); 

thr = 2.326; % one-tailed (p<0.01, 2.326); (p<0.05, 1.654);
sz = 16;

figure('Color','white')
namestr = ["Spikes(Bayesian)","Spikes(OLE)","FPA_{uhf}(OLE)"];
subplot(1,3,1)
is1 = 1; is2 =2;
x = wsZs(:,is1) > thr & wsZs(:,is2) > thr;
y = wsZs(:,is1) > thr & wsZs(:,is2) < thr;
z = wsZs(:,is1) < thr & wsZs(:,is2) > thr;

scatter(wsZs(:,is1),wsZs(:,is2),sz,'k','filled','MarkerEdgeColor','none'); hold on;

scatter(wsZs(x,is1),wsZs(x,is2),sz,'g','filled','MarkerEdgeColor','none');
scatter(wsZs(y,is1),wsZs(y,is2),sz,'b','filled','MarkerEdgeColor','none');
scatter(wsZs(z,is1),wsZs(z,is2),sz,'r','filled','MarkerEdgeColor','none');

plot([thr,thr],[-6,6]);
plot([-6,6],[thr,thr]);

xlabel(namestr(is1));ylabel(namestr(is2));
xlim([-2,6]);ylim([-2,6]); box on;

subplot(1,3,2)
is1 = 1; is2 = 3;
x = wsZs(:,is1) > thr & wsZs(:,is2) > thr;
y = wsZs(:,is1) > thr & wsZs(:,is2) < thr;
z = wsZs(:,is1) < thr & wsZs(:,is2) > thr;

scatter(wsZs(:,is1),wsZs(:,is2),sz,'k','filled','MarkerEdgeColor','none'); hold on;

scatter(wsZs(x,is1),wsZs(x,is2),sz,'g','filled','MarkerEdgeColor','none');
scatter(wsZs(y,is1),wsZs(y,is2),sz,'b','filled','MarkerEdgeColor','none');
scatter(wsZs(z,is1),wsZs(z,is2),sz,'r','filled','MarkerEdgeColor','none');

plot([thr,thr],[-6,6]);
plot([-6,6],[thr,thr]);

xlabel(namestr(is1));ylabel(namestr(is2));
xlim([-2,6]);ylim([-2,6]); box on;

subplot(1,3,3)
is1 = 2; is2 = 3;
x = wsZs(:,is1) > thr & wsZs(:,is2) > thr;
y = wsZs(:,is1) > thr & wsZs(:,is2) < thr;
z = wsZs(:,is1) < thr & wsZs(:,is2) > thr;

scatter(wsZs(:,is1),wsZs(:,is2),sz,'k','filled','MarkerEdgeColor','none'); hold on;

scatter(wsZs(x,is1),wsZs(x,is2),sz,'g','filled','MarkerEdgeColor','none');
scatter(wsZs(y,is1),wsZs(y,is2),sz,'b','filled','MarkerEdgeColor','none');
scatter(wsZs(z,is1),wsZs(z,is2),sz,'r','filled','MarkerEdgeColor','none');

plot([thr,thr],[-6,6]);
plot([-6,6],[thr,thr]);

xlabel(namestr(is1));ylabel(namestr(is2));
xlim([-2,6]);ylim([-2,6]); box on;

savefig(gcf,'Figure_supplement_5b','compact');

