clear; clc;
warning('off');

%% plot
load('Figure3ab_results_data.mat')
figure('Color','white');
median_error_FPAuhf = median_error_FPAuhf(1:5,1:5);
subplot(1,3,1)
[m, n] = size(median_error_FPAuhf);
imagesc(median_error_FPAuhf); colormap('jet'); clim([0,25]); hold on;
title('FPA_{uhf}');set(gca, 'YDir', 'reverse');
xticks([1:m]); xticklabels(char(64+[1:m]'));
yticks([1:n]); yticklabels(char(64+[1:n]'));
xlabel('Training Session'); ylabel('Testing Session');

textStrings = num2str(median_error_FPAuhf(:), '%0.1f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:m);  % Create x and y coordinates for the strings
idxc = x == y;
hStrings = text(x(~idxc), y(~idxc), textStrings(~idxc), ...  % Plot the strings
                'HorizontalAlignment', 'center', 'Color', 'k');
hStrings = text(x(idxc), y(idxc), textStrings(idxc), ...  % Plot the strings
                'HorizontalAlignment', 'center', 'Color', 'w');

subplot(1,3,2)
median_error_dLFPTheta = median_error_dLFPTheta(1:5,1:5);
imagesc(median_error_dLFPTheta); colormap('jet'); clim([0,20]);
title('dLFP_{\theta}');
xticks([1:m]); xticklabels(char(64+[1:m]'));
yticks([1:n]); yticklabels(char(64+[1:n]'));
xlabel('Training Session'); ylabel('Testing Session');
set(gca, 'YDir', 'reverse');
textStrings = num2str(median_error_dLFPTheta(:), '%0.1f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:m);  % Create x and y coordinates for the strings
idxc = x == y;
hStrings = text(x(~idxc), y(~idxc), textStrings(~idxc), ...  % Plot the strings
                'HorizontalAlignment', 'center', 'Color', 'k');
hStrings = text(x(idxc), y(idxc), textStrings(idxc), ...  % Plot the strings
                'HorizontalAlignment', 'center', 'Color', 'w');
            
subplot(1,3,3)
median_error_joint = median_error_joint(1:5,1:5);
imagesc(median_error_joint); colormap('jet'); clim([0,20]);
title('FPA_{uhf}+dLFP_{\Theta}');
xticks([1:m]); xticklabels(char(64+[1:m]'));
yticks([1:n]); yticklabels(char(64+[1:n]'));
xlabel('Training Session'); ylabel('Testing Session');
set(gca, 'YDir', 'reverse');
textStrings = num2str(median_error_joint(:), '%0.1f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:m);  % Create x and y coordinates for the strings
idxc = x == y;
hStrings = text(x(~idxc), y(~idxc), textStrings(~idxc), ...  % Plot the strings
                'HorizontalAlignment', 'center', 'Color', 'k');
hStrings = text(x(idxc), y(idxc), textStrings(idxc), ...  % Plot the strings
                'HorizontalAlignment', 'center', 'Color', 'w');
            
savefig(gcf,'Figure3a.fig','compact');

%%
figure('Color','white');
FPA = median_error_FPAuhf;
atheta = median_error_dLFPTheta;
joint = median_error_joint;
diagre = zeros(0);
for dd = -4:4
    deagre(dd+5,1) = mean(diag(FPA,dd));
    deagre(dd+5,2) = mean(diag(atheta,dd));
    deagre(dd+5,3) = mean(diag(joint,dd));
end

plot(deagre,'LineWidth',2); xticks([1,9]);
xticks([1,3,5,7,9]); xticklabels([-4,-2,0,2,4]); ylim([0,20]); yticks([0:5:20]);
legend({'FPA_{uhf}','dLFP_{\theta}','FPA_{uhf}+dLFP_{\theta}'},'Location','north');
ylabel('Decoding median error (cm)');
xlabel('k-th diagonal');
savefig(gcf,'Figure3b.fig','compact');
