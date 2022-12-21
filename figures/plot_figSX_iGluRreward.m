%% INPUTS
nAn = length(compare(1).s);
lickWithin = 0.25; %CHANGE, lick within this window
winRew = [-1 2]; % CHANGE, window for aligning signal to rewarded trials
winAcc = [-1 1]; % CHANGE, window for aligning signal to acceleration
NumStd = 2; % CHANGE, for immobility trough/peak analysis

for z = 1:length(compare)
    beh = compare(z).s;
    for x = 1:length(beh)
        Fs = beh(x).Fs;
        rew = beh(x).reward./Fs; % reward delivery times, in seconds
        lick = beh(x).lick./Fs; % lick times, in seconds
        [rewYes, ~, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
        %rewYes: indices of rewarded trials (animal licked within specified window)
        ev = rew(rewYes); % delivery times for rewarded trials, in seconds
        ev(ev < 1200) = []; ev(ev > 2400) = []; % post-infusion window
        for y = 1:length(beh(x).FP)
            sig = beh(x).FP{y}; % signal that will be aligned to event times
            sig = sig - nanmean(sig); % subtract mean of trace to center on zero
            [mat, time] = getSTA(sig, ev, Fs, [winRew(1), winRew(end)]);
            compare(z).a_rew{x,y} = mat; % save into structure
        end
    end
end

lag = []; val = []; % initialize output
y = 1; % analyzing ACh reward response
r = [200 700]; % range for ACh trough, in milliseconds
r = (r/(1000/Fs) + find(time == 0)); 
for z = 1:length(compare)
    for x = 1:size(compare(z).a_rew,1)
        pull = compare(z).a_rew{x,y}; % extract matrix of reward-aligned signals during post-infusion window into workspace
        pull = pull(r(1):r(2),:); % extract only segment of time specified as window for reward-related trough
        [a,b] = min(pull); % find local minima within range
        c = time(b + r(1) - 1); % convert index to seconds
        c(isnan(a)) = nan; 
        val(x,z) = nanmean(a); % save amplitude of minima
        lag(x,z) = nanmean(c).*1000; % save lag at minima, in milliseconds
    end
end

%%
fig = figure; hold on
plot(val','--.k','MarkerSize',20);
errorbar([0.75 2.25],nanmean(val),SEM(val,1),'.r','MarkerSize',20);
ylabel('ACh trough amp'); ylim([-5 0]); yticks([-8:2:0]);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'saline','iGluR'});
[~,p] = ttest2(val(:,1),val(:,2));
title(sprintf('Saline v iGluR antag: p = %1.2f (N = %d)',p,nAn));
axis square
movegui(gcf,'center');