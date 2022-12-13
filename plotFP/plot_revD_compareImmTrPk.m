%% LOAD DATA
% Comparing immobility ACh pause peak statistics
% for 2 different sets of recordings

<<<<<<< Updated upstream
beh1 = extractBeh; % Load PREVIOUS reward recordings
beh2 = extractBeh; % Load NEW reward recordings
compare = struct;
compare(1).s = beh1; compare(1).lbl = 'pre';
compare(2).s = beh2; compare(1).lbl = 'post';
=======
% beh1 = extractBeh; % Load PREVIOUS reward recordings
% beh2 = extractBeh; % Load NEW reward recordings
% compare = struct;
% compare(1).s = beh1; compare(1).lbl = 'pre';
% compare(2).s = beh2; compare(2).lbl = 'post';
>>>>>>> Stashed changes

%% INPUTS
nAn = length(compare(1).s);
% lickWithin = 0.25; %CHANGE, lick within this window
% winRew = [-1 2]; % CHANGE, window for aligning signal to rewarded trials
% winAcc = [-1 1]; % CHANGE, window for aligning signal to acceleration
NumStd = 2; % CHANGE, for immobility trough/peak analysis

%% Immobility pause peak
% NumStd = 2;
fig = figure; fig.Position([3 4]) = [1375 800];
lbl1 = {'trough','peak'};
lbl2 = {'frequency','duration','amplitude'};

for y = 1:length(compare)
    beh = compare(y).s; % Extract behavior
    [amp, dur, freq, thres] = getImmPausePeak(beh); % Immobility pause peak statistics
    Fs = beh(1).Fs;
    freq(isnan(freq)) = nan;
    dur(isnan(dur)) = nan; dur = dur.*(1000/Fs); % Adjust from samples to ms
    amp(isnan(amp)) = nan; amp = abs(amp); % Adjust to be absolute amplitude
    plotme = cell(1,3); plotme{1} = freq; plotme{2} = dur; plotme{3} = amp;
    cannula(y).out = plotme;

    nAn = length(beh);
    for jj = 1:2 % Iterate over pause and peaks
        for ii = 1:length(plotme)
            a = plotme{ii}(:,jj);
            subplot(2,3,ii + (jj-1)*3); hold on
            jit = []; 
            j1 = y-0.2; j2 = y+0.2; jit = j1 + (j2-j1).*rand(nAn,1); % jitter for pause
            plot(jit,a,'.k','MarkerSize',20); % plot raw values per animal
            errorbar([1,2],nanmean(a),SEM(a,1),'.r','MarkerSize',20); % plot average with error bars across animals
            xlim([0.5 2.5]); xticks([1 2]); xticklabels({compare.lbl});
            title(sprintf('%s-%s',lbl1{jj},lbl2{ii}));
            axis square;
        end
    end
end
fprintf('Done. \n'); 

lims1 = [0 1; 0 400; 0 6];
ticks1 = {[0:0.5:1], [0:100:500], [0:2:6]};
for jj = 1:2
    for ii = 1:3
        subplot(2,3,ii + (jj-1)*3);
        a = [cannula(1).plotme{ii}(:,jj), cannula(2).plotme{ii}(:,jj)];
        ylabel(sprintf('%s',lbl2{ii})); ylim(lims1(ii,:)); yticks(ticks1{ii});
        [~,p] = ttest(a(:,1),a(:,2));
        title(sprintf('%s-%s (p = %1.3f)',lbl1{jj},lbl2{ii}),p);
    end
end
movegui(gcf,'center');
