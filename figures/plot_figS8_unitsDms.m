%% LOAD
load('C:\Users\Anya\Desktop\FP_LOCAL\fig1\krok_figS8_unitsDMS.mat')

%% S8 B-E: comparing unit properties between DLS and DMS pCINs
fig = figure; fig.Position(3) = 1000;
lbl = {'pSPN','pCIN','other'}; % unit groups
clr = {'k','b'};
lblPlot = {'Firing rate (sp/s)','Coefficient of variation','Phasic activity index','Spike duration (ms)'};
lims = [0 15; 0 6; 0 0.3; 0 4];
for x = 1:length(units)
    vals = {[units(x).fr],[units(x).CV],[units(x).phasicIndex],[units(x).spkDur]};
    nUnits = length(vals{1});
    j1 = -0.3+x; j2 = 0.3+x; jit = j1 + (j2-j1).*rand(nUnits,1); % jitter
    for y = 1:4
        subplot(1,4,y); hold on
        plot(jit, vals{y}, '.', 'MarkerSize', 20, 'Color', clr{x});
        errorbar(x, nanmean(vals{y}), SEM(vals{y},2), '.r', 'MarkerSize', 20);
        ylabel(lblPlot{y}); ylim(lims(y,:));
        xlim([0.5 2.5]); xticks([1 2]);
        set(gca,'TickDir','out');
    end
end
y = 1; [~,p] = ttest2(units(1).fr, units(2).fr);
subplot(1,4,y); xticklabels({units.lbl}); title(sprintf('p = %1.3f',p));
y = 2; [~,p] = ttest2(units(1).CV, units(2).CV);
subplot(1,4,y); xticklabels({units.lbl}); title(sprintf('p = %1.3f',p));
y = 3; [~,p] = ttest2(units(1).phasicIndex, units(2).phasicIndex);
subplot(1,4,y); xticklabels({units.lbl}); title(sprintf('p = %1.3f',p));
y = 4; [~,p] = ttest2(units(1).spkDur, units(2).spkDur);
subplot(1,4,y); xticklabels({units.lbl}); title(sprintf('p = %1.3f',p));
movegui(gcf,'center');

%% S8 F-I: DMS pCIN synchrony
x = 1; % immobility
time = outCcgCcg(x).time;
nPairs = size(outCcgCcg(x).delta,2);

fig = figure;
fig.Position([3 4]) = [1000 800];
subplot(2,2,1); hold on
[~, ii] = sort(outCcg(x).delta(time == 0,:)); % sort in ascending order value at lag = 0
t = time(find(time == -1):find(time == 1));
h = imagesc(t, [1:nPairs], 1+outCcgCcg(x).delta(:,ii)', [0.8 1.5]);
colorbar; colormap(jet(256));
xlabel('Time to spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('pCIN pair');
h(x) = h.Parent; h(x).CLim = [0.8 1.5];
title(sprintf('CCG %s (n = %d)',outCcgCcg(x).lbl,nPairs));
axis square; set(gca,'TickDir','outCcg');

subplot(2,2,2); hold on
sm = 1;
shadederrbar(time, 1+movmean(nanmean(outCcg(x).delta50,2),sm), movmean(nanmean(outCcg(x).delta95,2),sm), 'k');
shadederrbar(time, 1+movmean(nanmean(outCcg(x).delta,2),sm), movmean(SEM(outCcg(x).delta,2),sm), 'b'); 
xlabel('Time to spike (s)'); xlim([-1 1]); 
ylabel('Firing rate (norm.)'); ylim([0.85 1.3]); yticks([0:0.1:1.4]);
legend({'shuff','immobility'});
title(sprintf('CCG mean, imm max = %1.2f',max(1+nanmean(outCcg(x).delta,2))));
axis square; set(gca,'TickDir','outCcg');

subplot(2,2,3); hold on
ds = 2;
a = 100*sum(outCcg(x).above95,2)/size(outCcg(x).above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(1:ds:end);
b = -100*sum(outCcg(x).below5,2)/size(outCcg(x).below5,2); b = b(1:ds:end);
bar(time(1:ds:end), a,'FaceColor','b','FaceAlpha',0.5);
bar(time(1:ds:end), b,'FaceColor','b','FaceAlpha',0.5);
xlabel('Time to spike (s)'); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('% of CCGs'); ylim([-50 100]); xlim([-1 1]);
legend({'immobility','immobility'});
title(sprintf('max = %1.1f @ %d ms',prop_m, round(1000*prop_t)))
axis square; set(gca,'TickDir','outCcg');

subplot(2,2,4); hold on
sm = 1;
x = 1; shadederrbar(time, 1+movmean(nanmean(outCcg(x).delta,2),sm), movmean(SEM(outCcg(x).delta,2),sm), 'k'); 
x = 2; shadederrbar(time, 1+movmean(nanmean(outCcg(x).delta,2),sm), movmean(SEM(outCcg(x).delta,2),sm), 'b'); 
xlabel('Time to spike (s)'); xlim([-1 1]); 
ylabel('Firing rate (norm.)'); ylim([0.9 1.4]); yticks([0:0.1:1.4]);
legend({outCcg.lbl});
title(sprintf('CCG mean, loc max = %1.2f',max(1+nanmean(outCcg(2).delta,2))));
axis square; set(gca,'TickDir','outCcg');

movegui(gcf,'center');

%% S8 K-L: Firing rate of DMS pCINs to ACh-DLS peaks
fig = figure; fig.Position(3) = 1000;
time = outACh(x).time;
nUnits = size(outACh(x).delta,2);
x = 1; % aligning to peaks

subplot(1,2,1); hold on
delta = out(x).delta - nanmean(outACh(x).delta([1:find(time == -0.51)],:)); % subtract baseline
shuff = out(x).delta50 - nanmean(outACh(x).delta50([1:find(time == -0.51)],:)); % subtract baseline
plot([0 0],[0.8 1.2],'--','Color',[0.5 0.5 0.5]);
shadederrbar(time, 1+movmean(nanmean(shuff,2),sm), movmean(nanmean(out(x).delta95,2),sm), 'k');
shadederrbar(time, 1+movmean(nanmean(delta,2),1), movmean(SEM(delta,2),1), 'b'); 
ylabel('pCIN DMS Firing rate (norm.)'); ylim([0.8 1.2]); yticks([0:0.1:2]);
xlabel(sprintf('Time to ACh-DLS %s (s)',outACh(x).lbl)); xlim([-1 1]); xticks([-1:0.5:1]);
title(sprintf('max = %1.3f',max(1+nanmean(outACh(x).delta,2))));
axis('square'); set(gca,'TickDir','out');

subplot(1,2,2); hold on
ds = 2;
a = 100*sum(outACh(x).above95,2)/size(outACh(x).above95,2); 
[prop_m,prop_t] = max(a); prop_t = time(prop_t); 
a = a(2:ds:end);
b = -100*sum(outACh(x).below5,2)/size(outACh(x).below5,2); b = b(1:ds:end);
bar(time(2:ds:end), a,'FaceColor','b','FaceAlpha',0.5,'EdgeAlpha',0.1);
bar(time(2:ds:end), b,'FaceColor','b','FaceAlpha',0.5,'EdgeAlpha',0.1);
xlabel(sprintf('Time to ACh-DLS %s (s)',outACh(x).lbl)); xlim([-1 1]); xticks([-1:0.5:1]);
ylabel('% of units'); ylim([-100 100]); yticks([-100:50:100]);
title(sprintf('max %1.1f @ %d ms',prop_m, round(1000*prop_t)))
axis('square'); set(gca,'TickDir','out');

movegui(gcf,'center');
%%
% x = 1;
% outCcg(x).delta = ccgDelta_rest;
% outCcg(x).delta50 = ccgDelta_50rest;
% outCcg(x).delta95 = ccgDelta_95rest;
% outCcg(x).above95 = above95;
% outCcg(x).below5 = below5;
% outCcg(x).time = time(:);
% 
% x = 2;
% outCcg(x).delta = ccgDelta_mvmt;
% outCcg(x).delta50 = ccgDelta_50mvmt;
% outCcg(x).delta95 = ccgDelta_95mvmt;
% outCcg(x).above95 = above95;
% outCcg(x).below5 = below5;
% outCcg(x).time = time(:);
