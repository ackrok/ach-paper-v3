load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figSX_naiveVwater_corr.mat');
% AK_corrFP

%%
fig = figure; fig.Position(3) = 1000;

subplot(1,2,1); hold on
plot(lags/Fs, corr_naive(:,5), 'k');
plot(lags/Fs, corr_water(:,5), 'b');
xlim([-2 2]); xticks([-2:1:2]); xlabel('Lag to DA (s)');
ylim([-0.5 0.1]); yticks([-1:0.25:0]); ylabel('coefficient');
legend({'naive','water-dep'},'Location','southeast');
title('matched ACh/DA cross-corr'); axis square

subplot(1,2,2); hold on
a = [min_val_naive, min_val_water];
plot(a', '.-k', 'MarkerSize', 20);
errorbar([0.75 2.25], nanmean(a), SEM(a,1), '.b', 'MarkerSize', 20);
ylim([-1 0]); yticks([-1:0.25:0]); ylabel('coefficient');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'naive','water-dep'});
[~,p] = ttest(a(:,1),a(:,2));
title(sprintf('p = %1.2f naive vs waterdep (n = %d)',p,size(a,1))); axis square