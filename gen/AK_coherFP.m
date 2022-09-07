function [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh)
%Compute coherogram between ACh and DA photometry signals acquired during
%dual color photometry recordings using Buzsaki function 'MTCoherogram'
%
%   [coher_achda, phase_achda, t, f] = AK_coherFP(beh)
%   [coher_achda, phase_achda, t, f, coher_shuff, phase_shuff] = AK_coherFP(beh)
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
%
%   Author: Anya Krok, March 2022

%% INPUTS
aa = [2 1]; % Photometry signal to use as reference is the one listed first
% e.g. if y = [1 2] then the signal in beh(x).FP{1} will be used as
% reference signal, while if y = [2 1] then the signal in beh(x).FP{2} will
% be used as reference signal instead
% Default from Krok 2022 is to use y = [2 1] so that rDA1m photometry
% signal is the reference signal
nStates = 3; % Number of behavioral states

%%
mat = struct;
h = waitbar(0,'coherogram');
for x = 1:length(beh)  % iterate over all recordings
  
    %% extract signals
    mat(x).rec = beh(x).rec; % load recording name
    fp_mat = [];
    fp_mat(:,1) = beh(x).FP{aa(1)}; Fs = beh(x).Fs; % extract photometry signal from structure, which will be used as reference
    fp_mat(:,1) = fp_mat(:,1) - nanmean(fp_mat(:,1)); % subtract baseline (mean of entire photometry signal) from fp
    fp_mat(:,2) = beh(x).FP{aa(2)}; % extract photometry signal from structure
    fp_mat(:,2) = fp_mat(:,2) - nanmean(fp_mat(:,2)); % subtract baseline (mean of entire photometry signal) from fp
    
    %% extract indices for behavioral states
    if isfield(beh,'reward')
        idx_rew = extractEventST([1:length(fp_mat)]', floor(beh(x).reward), floor(beh(x).reward)+(Fs*2), 1); % identify sample during reward
    else; idx_rew = []; end
    idx_mov = extractEventST([1:length(fp_mat)]', beh(x).on, beh(x).off, 1); % identify sample during locomotion
    idx_mov_nonRew = idx_mov(~ismember(idx_mov, idx_rew)); % exclude reward, include locomotion
    idx_imm = extractEventST([1:length(fp_mat)]', beh(x).onRest, beh(x).offRest, 1); % identify sample during locomotion
    idx_imm_nonRew = idx_imm(~ismember(idx_imm, idx_rew)); % exclude reward, include rest
    idx_c = cell(nStates,1); idx_c{1} = idx_imm_nonRew; idx_c{2} = idx_mov_nonRew; idx_c{3} = idx_rew; % index into cell array for ease of iteration
    %%
    for z = 1:nStates % iterate over behavioral states
        if ~isempty(idx_c{z})
            sig = fp_mat(idx_c{z},:); % extract indexes samples 
            [coher,ph,t,f] = bz_MTCoherogram(sig(:,1),sig(:,2),'frequency',Fs,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
            mat(x).coher(:,z) = nanmean(coher,2); % collapse time dimension
            mat(x).phase(:,z) = nanmean(ph,2); % collapse time dimension
        end
    end
    %%
    tmp_coher = []; tmp_phase = []; 
    fp_shuff = fp_mat(:,2);
    for s = 1:50 % repeat shuffle N times
        % fp_shuff = circshift(fp_shuff, Fs);
        fp_shuff = fp_shuff(randperm(length(fp_shuff)));
        [idx_c,p] = bz_MTCoherogram(fp_mat(:,1),fp_shuff,'frequency',Fs,'range',[0 10],'window',10,'overlap',5,'step',5,'tapers',[3 5],'pad',0);
        tmp_coher(:,s) = nanmean(idx_c,2); % collapse time dimension
        tmp_phase(:,s) = nanmean(p,2); % collapse time dimension
    end
    mat(x).coher_shuff = prctile(tmp_coher, [5 50 95], 2); % 5th, 50th, 95th percentiles
    mat(x).phase_shuff = prctile(tmp_phase, [5 50 95], 2); % 5th, 50th, 95th percentiles
    
    waitbar(x/length(beh),h);
end
close(h); fprintf('Coherogram Analysis Done !\n');

%% EXTRACT COHERENCE and PHASE DURING EACH BEHAVIORAL STATE
coher_state = cell(nStates,1); 
phase_state = cell(nStates,1);
for x = 1:length(mat) % iterate over all recordings
    for z = 1:nStates % iterate over behavioral states
        if z == 3 && size(mat(x).coher,2) < 3 % if no coherence or phase for reward
            coher_state{z}(:,x) = nan(103,1); 
            phase_state{z}(:,x) = nan(103,1);
        else
        coher_state{z}(:,x) = mat(x).coher(:,z); % concatenate
        phase_state{z}(:,x) = mat(x).phase(:,z);
        end
    end
end

%% AVERAGE ACROSS ALL RECORDINGS FOR ONE ANIMAL SUCH THAT N = X mice
rec = {}; for x = 1:length(beh); rec{x} = strtok(beh(x).rec,'-'); end
uni = unique(rec); 
nAn = length(uni); % number of unique animal IDs

coher_achda = cell(3,1); phase_achda = cell(3,1); % initialize output
coher_shuff = cell(3,1); phase_shuff = cell(3,1); % initialize output

for x = 1:nAn
    ii = find(strcmp(rec,uni{x})); % index of matching recordings for this animal
    coh_tmp = [mat(ii).coher_shuff];
    ph_tmp = [mat(ii).phase_shuff];
    for z = 1:nStates
        coher_achda{z}(:,x) = nanmean(coher_state{z}(:,ii),2); % average recordings for each animal, for each behavioral state
        phase_achda{z}(:,x) = nanmean(phase_state{z}(:,ii),2); % average recordings for each animal, for each behavioral state
        coher_shuff{z}(:,x) = nanmean(coh_tmp(:,[z:3:size(coh_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
        phase_shuff{z}(:,x) = nanmean(ph_tmp(:,[z:3:size(ph_tmp,2)]),2); % average recordings for each animal, for each shuff percentile
    end
end
