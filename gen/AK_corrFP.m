function [corr_achda, lags, shuff_achda] = AK_corrFP(beh)
%Cross-correlation between ACh and DA photometry signals acquired during
%dual color photometry recordings
%
%   [corr_achda, lags] = AK_corrFP(beh)
%   [corr_achda, lags, shuff_achda, min_val, min_lag] = AK_corrFP(beh)
%
%   Description: This function is for running cross-correlation analysis on
%   two continuous photometry signals (ACh, DA) after isolating time
%   periods when animal is immobile (no movement, no reward) using the
%   MATLAB function 'xcorr'
%
%   INPUTS
%   'beh' - structure with photometry and behavioral data for multiple
%   recordings, should include beh(x).rec as [an,'-',day]
%   
%   OUTPUTS
%   'corr_achda' - cross-correlation between ACh and DA signals, normalized
%   using coeff scaling to approximate Pearson's coefficient, r.
%   'lags' - vector of lag indices from +/- 5 seconds, in samples
%   'shuff_achda' - 5th, 50th, 95th percentiles for cross-correlation run
%   on shuffled photometry signals, across behavioral states
%       shuff_achda{1,1} is 5th percentile for immobility, and so on
%       shuff_achda{3,2} is 50th percentile for reward, and so on
%       shuff_achda{2,3} is 95th percentile for locomotion, and so on
%
%   Author: Anya Krok, December 2021

%% INPUTS
y = [2 1]; % Photometry signal to use as reference is the one listed first
% e.g. if y = [1 2] then the signal in beh(x).FP{1} will be used as
% reference signal, while if y = [2 1] then the signal in beh(x).FP{2} will
% be used as reference signal instead
% Default from Krok 2022 is to use y = [2 1] so that rDA1m photometry
% signal is the reference signal

%% OUTPUTS
% mat = struct; % temporary output structure
corr_cell = cell(3,4); % temporary output cell array
for a = 1:3; for b = 1:4; corr_cell{a,b} = nan(501,length(beh)); end; end % fill cell array with nan's

%% RUN ANALYSIS ON ALL RECORDINGS
h = waitbar(0, 'cross correlation');
for x = 1:length(beh) % iterate over all recordings
    
    %% extract signals
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{y(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
    fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
    fp_mat(:,2) = beh(x).FP{y(2)}; % extract photometry signal from structure
    fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp

    %% extract indices for behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+Fs, 1); % identify sample during reward window: [0 +1]seconds relative to reward delivery for all trials (rewarded or not)
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_cell = cell(3,1); idx_cell{1} = idx_imm_nonRew; idx_cell{2} = idx_mov_nonRew; idx_cell{3} = idx_rew; % index into cell array for ease of iteration
    
    %%
    % mat(x).rec = beh(x).rec;
    
    fp_sub = [];
    for z = 1:3
        if length(idx_cell{z})< 2; continue; end
        fp_sub = fp_mat(idx_cell{z},:); % signal
        % fp_sub = normalize(fp_mat,1,'range'); % normalize [0 1]
        % fp_sub = fp_filt(idx_cell{z},:); % bandpass filtered signal
        
        % [xcf, lags, bounds] = crosscorr(fp_sub(:,1), fp_sub(:,2),'NumLags',100,'NumSTD',3);
        % [shuff,~,~] = crosscorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2),'NumLags',100,'NumSTD',3);
        [corr_tmp, lags] = xcorr(fp_sub(:,1), fp_sub(:,2), 5*Fs, 'coeff'); % cross-correlation
        
        fp_sub_new = fp_sub(:,2);
        tmp_shuff = []; 
        for s = 1:50
            fp_sub_new = circshift(fp_sub_new, Fs);
            % tmp_shuff(:,s) = xcorr(fp_sub(randperm(size(fp_sub,1)),1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
            % tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub(randperm(size(fp_sub,2)),2), 10*Fs, 'coeff');
            tmp_shuff(:,s) = xcorr(fp_sub(:,1), fp_sub_new, 5*Fs, 'coeff');
        end
        
        % mat(x).corr{z} = corr_tmp;
        % mat(x).shuff{z} = prctile(tmp_shuff, [5 50 95], 2);
        corr_cell{z,1}(:,x) = corr_tmp;       % cross-correlation
        corr_cell{z,2}(:,x) = prctile(tmp_shuff, 5, 2); % shuffle 5th percentile
        corr_cell{z,3}(:,x) = prctile(tmp_shuff, 50, 2); % shuffle 50th percentile
        corr_cell{z,4}(:,x) = prctile(tmp_shuff, 95, 2); % shuffle 95th percentile
    end
    
%%
    waitbar(x/length(beh),h);
end
close(h);

%% AVERAGE ACROSS ALL RECORDINGS FOR ONE ANIMAL SUCH THAT N = X mice
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); 
nAn = length(uni); % number of unique animal IDs

corr_achda = cell(3,1); % initialize output
shuff_achda = cell(3,3); % initialize output 

for x = 1:nAn
    idx = strcmp(rec,uni{x}); % match animal ID to recordings
    for z = 1:3 % iterate over behavioral states
        corr_adj = corr_cell{z,1}; % extract cross-correlation output for this behavioral state
        corr_adj = corr_adj - nanmean(corr_adj([1:find(lags == -2)],:)); % adjust such that baseline outside of +/- 2s is at zero
        corr_achda{z,1}(:,x) = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
        for b = 2:4 % iterate over shuffle percentiles
            corr_adj = corr_cell{z,b};
            shuff_achda{z,b-1} = nanmean(corr_adj(:,idx),2); % average across all recordings for this animal
        end
    end
end

end