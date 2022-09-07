behLoaded = menu('beh loaded into workspace already?','yes','no');

switch behLoaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        if ~iscell(fName); fName = {fName}; end

        %% Extract data
        if ~exist('beh','var'); beh = struct; end
        for y = 1:length(fName) 
            load(fullfile(fPath,fName{y})); 
            [an,b] = strtok(fName{y},'_'); day = strtok(b,'_');
            x = 1+length(beh);
            beh(x).rec = [an,'-',day]; 
            beh(x).site = 'DLS';
            beh(x).task = 'wheel';

            %% Photometry
            beh(x).Fs = data.gen.Fs; % Sampling frequency, in Hz
            beh(x).time = data.final.time; % Time vector
        %     beh(x).FP = data.final.FP;
        %     beh(x).nbFP = data.final.nbFP; % Photometry signal(s)
        %     beh(x).FPnames = data.final.FPnames; % Names of photometry signal(s)

            %% Movement
            if isfield(data.final,'vel') % If movement data exists
                beh(x).vel = data.final.vel; % Velocity signal
                beh(x).on = data.final.mov.onsets; beh(x).off = data.final.mov.offsets;                 % Movement onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).onRest = data.final.mov.onsetsRest; beh(x).offRest = data.final.mov.offsetsRest; % Rest onset/offset times in sampling freq (data.gen.Fs), NOT in seconds
            end

            %% Lick/Reward
            if isfield(data.acq,'rew')
                % data = processReward(data, data.gen.params);
                beh(x).task = 'reward';
                beh(x).reward = data.final.rew.onset;    % Reward delivery time in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).lick = data.final.lick.onset;        % Lick times in sampling freq (data.gen.Fs), NOT in seconds
                beh(x).lickVec = data.final.lick.trace;        % Lick trace
            end

            %%
            fprintf('Extracted from %s\n',fName{y});
        end
        if isempty(beh(1).Fs); beh(1) = []; end
end

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