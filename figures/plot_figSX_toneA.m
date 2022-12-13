% figure;
% clr = {'k','r','b','g','m'};
% for x = 1:4
%     subplot(2,3,x); hold on
%     for z = 1:4
%         plot(c1(z).s(x).time, c1(z).s(x).nbFP{2}, clr{z});
%     end
%     title(sprintf('%s',strtok(c1(1).s(x).rec,'-')));
% end
% %%
a = cell(1,2);
c1 = cannula; c1(3) = [];
winInf = [0 40]; Fs = 50;
for z = 1:length(c1)
    beh = c1(z).s;
    for x = 1:length(beh)
        if isempty(beh(x).nbFP); a{1}(x,z) = nan; a{2}(x,z) = nan; continue; end
        for y = 1:2
            demod = beh(x).nbFP{y};
            idxImm = extractEventST([1:length(demod)]', beh(x).onRest, beh(x).offRest,1); % Sample index during immobility
            idxImm(idxImm < winInf(1)*60*Fs | idxImm > winInf(2)*60*Fs) = []; % Remove sample indices that are outside of infusion window
            demodImm = demod(idxImm); % Restrict to immobility + infusion window
            demodBase = mode(demod(1:10000));
            a{y}(x,z) = min(demodImm)/demodBase;
        end
    end
end
a{1}([5:9], 4) = nan; a{2}([5:9], 4) = nan;
%% DA
global_da = a{2};
z = 4; global_da(:,z) = (global_da(:,z) - (global_da(:,z)))./(nanmean(global_da(:,1)) - nanmean(global_da(:,z)));
global_ach = a{1};
z = 3; global_ach(:,z) = (global_ach(:,z) - (global_ach(:,z)))./(nanmean(global_ach(:,1)) - nanmean(global_ach(:,z)));

%%
lbl = {'saline','iGluR anatg','mAChR antag','D1/2R antag'};

[~,p] = ttest(global_da(:,1), global_da(:,2));
[~,p(2)] = ttest(global_da(:,1), global_da(:,3));
[~,p(3)] = ttest(global_da(:,1), global_da(:,4));
fig = figure; fig.Position(3) = 1000;
subplot(1,2,2); hold on
plot(global_da','.m', 'MarkerSize', 20); 
errorbar([1.2:1:4.2], nanmean(global_da), SEM(global_da,1), '.k', 'MarkerSize', 20); 
xlim([0.5 4.5]); xticks([1:5]); xticklabels(lbl);
ylabel('Global DA signal'); yticks([0:0.25:1.25]); ylim([0 1.25]);
title(sprintf('DA: sal/iglur = %1.3f | sal/machr = %1.3f | sal/d1d2 = %1.5f',p(1),p(2),p(3)))
title(sprintf('DA: sal/iglur = %1.3f | sal/machr = %1.3f | sal/d1d2 = 1.1e-04',p(1),p(2)))

[~,p] = ttest(global_ach(:,1), global_ach(:,2));
[~,p(2)] = ttest(global_ach(:,1), global_ach(:,3));
[~,p(3)] = ttest(global_ach(:,1), global_ach(:,4));
subplot(1,2,1); hold on
plot(global_ach','.g', 'MarkerSize', 20); 
errorbar([1.2:1:4.2], nanmean(global_ach), SEM(global_ach,1), '.k', 'MarkerSize', 20); 
xlim([0.5 4.5]); xticks([1:5]); xticklabels(lbl);
ylabel('Global ACh signal'); yticks([0:0.25:1.25]); ylim([0 1.25]);
title(sprintf('ACh: sal/iglur = %1.3f | sal/machr = %1.3f | sal/d1d2 = %1.3f',p(1),p(2),p(3)))
title(sprintf('ACh: sal/iglur = %1.3f | sal/machr = 6.9e-06 | sal/d1d2 = %1.3f',p(1),p(3)))