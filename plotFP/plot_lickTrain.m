behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        beh = extractBeh(fPath, fName);
end
if ~exist(beh); error('No variable called beh exists'); end

%% Plot lick rate
lickWithin = 0.25; %CHANGE, lick within this window

figure; plm = floor(sqrt(length(beh))); pln = ceil(length(beh)/plm);
for x = 1:length(beh)
    Fs = beh(x).Fs;
    rew = beh(x).reward./Fs; % reward delivery times, in seconds
    lick = beh(x).lick./Fs; % lick times, in seconds
    [rewYes, rewNo, lickNew] = extractRewardedTrials(rew, lick, [0 lickWithin]);
    %rewYes: indices of rewarded trials (animal licked within specified window)
    %rewNo: indices of non-rewarded trials
    %lickNew: corrected vector of lick times, in seconds
    bin = 0.1; window = [-2 2];
    peth = getClusterPETH(lickNew, rew(rewYes), bin, window); % align lick to rewarded trials
    sp(x) = subplot(plm,pln,x); hold on
    plot([0 0],[0 1],'k');
    shadederrbar(peth.time, nanmean(peth.cts{1},2), SEM(peth.cts{1},2), 'm'); % plot lick rate
    xlabel('latency to reward (s)'); ylabel('lick rate');
    title(sprintf('%s: (%d/%d)',beh(x).rec,length(find(rewYes)),length(rew)),'Interpreter','none');
    %title of subplot is 
    % <recording name (#rewarded trials / total #trials)>
end
