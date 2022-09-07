%% Description
% The code below will use FFT to analyze the frequency components present
% in photometry signal with particular emphasis on delta frequency band
% found to be enriched in DLS ACh and DA photometry signals
%
% (1) Dialog window will prompt user to specify inputs. Options are in
%   parentheses!
% (2) If raw signals are not already loaded into workspace in a structure 
%   called rawS, user will then be prompted to select the .mat files that 
%   contain photometry data to be analyzed.
% (3) User will be prompted to select .mat file that contains normalized
%   FFT vectors for recordings of GFP and tdTomato.
% (4) FFT will be run on signal extracted during specified behavioral state,
%   then normalized between 0.1 and 100 Hz, and lastly the FFT power for 
%   GFP or tdTomato will be subtracted.
% (5) FFT output will be plotted if specified to do so.
%
% INPUTS
%   Specify in dialog window.
%   Background functioned needed to run this code: extractEventST
%   Data files needed: raw data files (.mat), GFP+tdTomato data file (.mat)
%
% OUTPUTS
%   'flog' - frequency vector (log-scale),  for plotting FFT output
%   'sub_mat' - matrix with FFT output for all recordings
%   'auc' - area under the curve in 0.5-4Hz frequency band
%
% Anya Krok, July 2022

%% INPUTS
behState = 'immobility'; % Options: 'immobility','locomotion','reward','full'
range_norm = [0.1 100]; % Range for normalization of FFT output, in Hz
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
fluorophore = 'green'; % Options: 'green','red' -- green: use GFP, red: use tdTomato for subtraction
y = 1; % Options: 1, 2 -- 1:ACh, 2:redDA
plotYes = 'yes'; % Options: 'yes','no' -- yes plot, no plot

answer = inputdlg({...
    'Behavioral state: (full, immobility, locomotion, reward)',...
    'Photometry signal: (1, 2)',...
    'Fluorophore: (green, red)',...
    'Plot FFT output? (yes, no)'},...
    'Inputs for FFT', [1 40],...
    {behState,num2str(y),fluorophore,plotYes});

behState = answer{1};
y = str2num(answer{2});
fluorophore = answer{3};
plotYes = answer{4};

%% LOAD RAW SIGNALS INTO WORKSPACE
loaded = menu('Already loaded raw data into workspace?','yes','no');
switch loaded
    case 2
        %% Select .mat files you want to add to summary data structu
        fPath = 'R:\tritsn01labspace\'; 
        [fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
        if ~iscell(fName); fName = {fName}; end

        %% Extract data
        rawS = struct;
        h = waitbar(0, 'new processing');
        for f = 1:length(fName)
            load(fullfile(fPath,fName{f})); % Load data file
            [an,b] = strtok(fName{f},'_'); day = strtok(b,'_'); % Parse file name
            x = 1 + length(rawS);
            rawS(x).rec = [an,'-',day]; 
            rawS(x).site = 'DLS';

            %% Pull parameters required for this analysis
            if isfield(data.gen,'params')
                params = data.gen.params; % Extract params structure
                dsRate = params.dsRate; 
                dsType = params.dsType; % General downsampling parameter
                rawFs = data.gen.acqFs; 
                Fs = data.gen.Fs;
            end
            %% Process photometry data
            fprintf('Extracting raw photometry data %s ... ',rawS(x).rec);
            rawS(x).FPnames = data.acq.FPnames;
            rawS(x).rawFP = data.acq.FP;
            fprintf('DONE.\n');
            rawS(x).rawFs = rawFs;
            if isfield(data,'final'); if isfield(data.final,'mov')
                rawS(x).on = data.final.mov.onsets.*dsRate;
                rawS(x).off = data.final.mov.offsets.*dsRate;
                rawS(x).onRest = data.final.mov.onsetsRest.*dsRate;
                rawS(x).offRest = data.final.mov.offsetsRest.*dsRate;
                end; end
            if isfield(data.final,'rew'); if isfield(data.final.rew,'onset')
                rawS(x).reward = data.final.rew.onset.*dsRate;
                end; end
            waitbar(f/length(fName),h);
        end
        close(h);
        if isempty(rawS(1).rawFs); rawS(1) = []; end
        
        %% extract during IMMOBILITY
        rmv = zeros(length(rawS),1);
        switch behState
            case 'immobility'
                fprintf('Extracting signal during %s (this will take a while!) ...',behState)
                h = waitbar(0,'Extracting signal during behavioral state');
                for x = 1:length(rawS)
                    nSampRaw = length(rawS(x).rawFP{y});
                    if isfield(rawS,'reward')
                        rewWindow = rawS(x).rawFs;
                        idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
                    else; idx_rew = [];
                    end
                    idx_imm = extractEventST([1:nSampRaw]', rawS(x).onRest, rawS(x).offRest, 1); % identify recording indices during immobility
                    if isempty(idx_imm); fprintf('%s - no immobility \n', rawS(x).rec); rmv(x) = 1; end
                    idx_imm = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward from immobility indices
                    rawS(x).fp_sub = rawS(x).rawFP{y}(idx_imm); % extract signal during immobility
                    waitbar(x/length(rawS),h);
                end; close(h);
            case 'locomotion'
                fprintf('Extracting signal during %s (this will take a while!) ...',behState)
                h = waitbar(0,'Extracting signal during behavioral state');
                for x = 1:length(rawS)
                    nSampRaw = length(rawS(x).rawFP{y});
                    if isfield(rawS,'reward')
                        rewWindow = rawS(x).rawFs;
                        idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
                    else; idx_rew = [];
                    end
                    idx_loc = extractEventST([1:nSampRaw]', rawS(x).on, rawS(x).off, 1); % identify recording indices during locomotion
                    if isempty(idx_loc); fprintf('%s - no locomotion \n', rawS(x).rec); rmv(x) = 1; end
                    idx_loc = idx_loc(~ismember(idx_loc, idx_rew)); % exclude reward from locomotion indices
                    rawS(x).fp_sub = rawS(x).rawFP{y}(idx_loc); % extract signal during locomotion
                    waitbar(x/length(rawS),h);
                end; close(h);
            case 'reward'
                fprintf('Extracting signal during %s (this will take a while!) ...',behState)
                h = waitbar(0,'Extracting signal during behavioral state');
                for x = 1:length(rawS)
                    nSampRaw = length(rawS(x).rawFP{y});
                    if isfield(rawS,'reward')
                        rewWindow = rawS(x).rawFs;
                        idx_rew = extractEventST([1:nSampRaw]', floor(rawS(x).reward), floor(rawS(x).reward)+rewWindow, 1); % identify recording indices during reward
                        rawS(x).fp_sub = rawS(x).rawFP{y}(idx_rew); % extract signal during reward
                    else
                        error('no reward');
                    end
                    waitbar(x/length(rawS),h);
                end; close(h);
            case 'full'
                for x = 1:length(rawS)
                    rawS(x).fp_sub = rawS(x).rawFP{y}; % full trace
                end
        end
        rawS(rmv == 1) = []; % remove recordings where no fp_sub extracted
        fprintf('DONE! \n');
end

%% FFT
p1_mat = [];
for x = 1:length(rawS)
    vec = [rawS(x).fp_sub]; 
    Fs = rawS(x).rawFs;
    
    needL = 2500*Fs;
    vec = repmat(vec,[ceil(needL/length(vec)) 1]);
    vec = vec(1:needL);
    T = 1/Fs;               % Sampling period
    L = length(vec);        % Length of signal
    vec(isnan(vec)) = [];
    fftACh = fft(vec);      % Discrete Fourier Transform of photometry signal
    P2 = abs(fftACh/L);     % Two-sided spectrum P2
    P1 = P2(1:L/2+1);       % Single-sided spectrum P1 based on P2 and even-valued signal length L
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;     % Frequency domain vector
    P1 = medfilt1(P1);      % Median filter initial FFT
    p1_mat(:,x) = [movmean(P1,500)];
end
fprintf('FFT done! (1/3) ... ');

%% Normalize FFT
tmp = [];
r = [find(f == range_norm(1)):find(f == range_norm(2))]; % Restrict to [0.01 100]
flog = log10(f(r));
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
for x = 1:size(p1_mat,2)
    a = log10(p1_mat(r,x));
    vec_norm = (a - a(end))./(a(1) - a(end)); 
    % vec_norm = normalize(log10(p1_mat(r,x)),'range'); % Normalize range from [0.01 100], scaling so range covers [0 1]
    tmp(:,x) = vec_norm;
end
norm = tmp;
fprintf('Normalization done! (2/3) ... ');

%% Subtract stable fluorophore (GFP/tdTomato) signal
[fName_flu,fPath_flu] = uigetfile([fPath,'*.mat'],'Select GFP+tdTomato for FFT file','MultiSelect','Off');
load(fullfile(fPath_flu,fName_flu));

sub_1 = [norm]; % SUBTRACT
switch fluorophore
    case 'green'; sub_2 = norm_gfp; % FFT ouput: GFP fluorescence signal, average over n = 3 mice
    % case 'red'; sub_2 = norm_tdt; % FFT ouput: tdTomato fluorescence signal, average over n = 3 mice
    case 'red'; sub_2 = norm_antd1d2; % FFT ouput: rDA fluorescence signal during infusion of DA receptor antagonist, average over n = 4 mice
end
sub_mat = []; for x = 1:size(sub_1,2); sub_mat(:,x) = sub_1(:,x) - nanmean(sub_2,2); end % Subtract avg FFT for mAChR antagonist
fprintf('Subtraction of %s fluorophore done! (3/3)\n', fluorophore);

%% AUC in specified frequency band
auc = [];
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
for x = 1:size(sub_mat,2)
    auc(x) = trapz(sub_mat(r_auc,x))/length(r_auc);
end

%% PLOT
switch plotYes
    case 'yes'
        ds = 50; % Downsample when plotting so figure is smaller in size
        nAn = length(rawS); % Number of recordings and/or animals
        fig = figure; fig.Position(3) = 1375;
        subplot(1,3,1); hold on
            plot(flog(1:ds:end), sub_mat((1:ds:end),:)); % Plot FFT for each recording
            plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
            plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
            plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
            xlabel('Frequency'); ylabel('Power (a.u.)');
            xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
            legend({rawS.rec});
            title(sprintf('FFT %s',behState)); axis square
        subplot(1,3,2); hold on
            shadederrbar(flog(1:ds:end), nanmean(sub_mat((1:ds:end),:),2), SEM(sub_mat((1:ds:end),:),2), 'r'); % Plot FFT averaged across all recordings
            plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
            plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
            plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
            xlabel('Frequency'); ylabel('Power (a.u.)');
            xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
            title(sprintf('FFT %s (n = %d)',behState,nAn)); axis square
        subplot(1,3,3); hold on
            j1 = 0.8; j2 = 1.2; jit = j1 + (j2-j1).*rand(nAn,1); % jitter
            plot(jit, auc, '.k', 'MarkerSize', 20); % Plot area under curve power data points for each recording
            errorbar(nanmean(auc), SEM(auc,2), '.r', 'MarkerSize', 20); % Plot area under curve power averaged
            xlim([0.5 1.5]); xticks([1]); xticklabels({'AUC'});
            ylabel('Power (a.u.)'); ylim([0 0.5]);
            title(sprintf('AUC %1.1f-%dHz (mu: %1.2f)',range_auc(1),range_auc(2),nanmean(auc))); axis square
        movegui(gcf,'center');
end