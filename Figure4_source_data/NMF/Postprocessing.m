clear; clc
warning('off');

load('pyhsmm_data_run.mat');
load('pyhsmm_data_run_result.mat');

figure('Color','white');
subplot(1,2,1)
N = histcounts2(rpos,fpa_train_inf',linspace(0,1,30),1:42);
[~,imax] = max(N,[],1);
[~,ord] = sort(imax);
imagesc(N(:,ord)'); colormap(flipud(bone));
xlabel('Real position'); ylabel('Stets infered from unsorted FPA_{uhf}');

subplot(1,2,2)
N = histcounts2(spk_train_inf,fpa_train_inf,1:max(spk_train_inf'),1:42);
[~,imax] = max(N,[],1);
[~,ord] = sort(imax);
imagesc(N(:,ord)'); colormap(flipud(bone));
xlabel('States infered from Sorted Spikes'); ylabel('States infered from unsorted FPA_{uhf}');
suptitle('Running')


clear; clc;
load('pyhsmm_data_ripple_result.mat');

figure('Color','white');
N = histcounts2(spk_train_inf,fpa_train_inf,1:max(spk_train_inf'),1:max(fpa_train_inf'));
[~,imax] = max(N,[],1);
[~,ord] = sort(imax);
imagesc(N(:,ord)'); colormap(flipud(bone));
xlabel('States infered from Sorted Spikes'); ylabel('States infered from unsorted FPA_{uhf}');
suptitle('Ripple events')

% spk_rates_mat = spk_rates_mat(1:spk_N_used_inf,:);
% [~,imax] = max(spk_rates_mat,[],1);
% [~,ord] = sort(imax);
% imagesc(spk_rates_mat(:,ord)'); colormap(flipud(bone));