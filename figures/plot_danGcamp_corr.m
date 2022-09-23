load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_gcamp-DLS-SNc');

%%
% figure;
% for x = 1:length(beh)
%     subplot(3,3,x); hold on
%     for y = 1:length(beh(x).FP)
%         plot(beh(x).FP{y});
%     end
% end

%% COMPUTE CROSS_CORRELATION
mat = struct;
for x = 1:length(beh); y = [2 1]; %CHANGE - which FP signal to run CCG or ACG on
    %% extract signals
    Fs = beh(x).Fs; 
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{y(1)}; % extract photometry signal from structure
    fp_mat(:,2) = beh(x).FP{y(2)};
    fp_mat = fp_mat - nanmean(fp_mat); % center traces
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    fp_sub = fp_mat(idx_imm,:); % extract signal during immobility
    [c,lags] = xcorr(fp_sub(:,1), fp_sub(:,2), 10*Fs, 'coeff'); % cross-correlation during immobility
    tmp_shuff = []; 
    for s = 1:50
        fp_shift = fp_sub(randperm(length(fp_sub)),2); % shuffle photometry values
        tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_shift, 10*Fs, 'coeff'); % cross-correlation of shuffled DLS signal
    end
    mat(x).rec = beh(x).rec;
    mat(x).ref = beh(x).FPnames{y(1)};
    mat(x).run = beh(x).FPnames{y(2)};
    mat(x).lags = lags./Fs;
    mat(x).c = c;
    mat(x).shuff = prctile(tmp_shuff, [2.5 50 97.5], 2);
end
fprintf('Cross-correlation GCaMP6f DLS/SNc DONE!\n');

%% EXTRACT MAXIMUM COEFFICIENT by mouse
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); nAn = length(uni);
mat_corr = []; % initialize output
for x = 1:nAn
    idx = find(strcmp(rec,uni{x})); % identify all recordings from this animal
    [~,idx2] = max(max([mat(idx).c])); % maximum amongst corr for all recordings
    mat_corr(:,x) = [mat(idx(idx2)).c]; % extract recording with highest correlation
end
[mm, ii] = max(mat_corr); % peak correlation coefficient
ii = lags(ii)/Fs*1000; % convert to milliseconds

%% PLOT
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
x = 4; % example recording
plot(lags/Fs, mat(x).c, 'm');
% shadederrbar(lags/Fs, mat(x).shuff(:,2), mat(x).shuff(:,1), 'k');
title('example corr imm DLS/SNc'); axis square
xlabel('Latency to SNc (s)');
ylabel('Coefficient'); ylim([-0.2 1]); yticks([-0.2:0.2:1]);

subplot(1,2,2); hold on
clr = hsv(nAn);
for x = 1:nAn
    plot(ii(x), mm(x), '.', 'MarkerSize', 20, 'Color', clr(x,:)); % plot individual data points
end
plot([0 0],[0 1],'--k'); % vertical line at lag = 0
errorbar(nanmean(ii), nanmean(mm), SEM(mm,2), '.k', 'MarkerSize', 20); % error bar with SEM
title(sprintf('Peak: %1.2f +/- %1.2f || Lag: %1.2f +/- %1.2f ms', nanmean(mm), nanstd(mm), nanmean(ii), nanstd(ii))); 
xlabel('Latency to SNc (ms)'); xlim([-400 400]); xticks([-500:100:500]);
ylabel('Coefficient'); ylim([0 1]); yticks([0:0.2:1]);
legend(uni);
axis square
movegui(gcf,'center');