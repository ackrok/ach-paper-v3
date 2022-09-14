%% Description
% The code below will use compare frequency bands enriched in specified
% photometry signals during specified behavioral state between two
% different intra-cranial infusions.
%
% INPUTS
%   'norm' - matrix with normalized FFT output for all recordings
%   'norm_cell' - cell array with normalized FFT output, with
%       each cell containing data for one infusion
%   'lbl_inf' - cell array with labels for infusions included in norm_cell
%   'lbl_fp' - cell array with FPname for photometry signal analyzed
%   'f' - frequency domain vector
%   'flog' - log-scale frequency domain vector
%
% OUTPUTS
%   'sub_comp' - cell array with FFT output for all recordings, 
%       each cell contains data for one infusion, subtracting
%       green or red fluorophore, as specified
%   'auc' - area under the curve in 0.5-4Hz frequency band
%       each column contains data for one infusion
%
% Anya Krok, July 2022

%%
if ~exist('rawS','var') && ~exist('norm','var')
    [norm, f, flog, p1_mat, rawS] = getCannula_fft();
end

%% INPUTS - MANUAL
% nComp = 2; % Number of comparisons to be made
% norm_comp = cell(1,nComp);
% norm_comp{1} = []; norm_comp{2} = []; % CHANGE
% lbl_inf = {'DLS','DMS'};
% lbl_fp = rawS(1).fp_lbl;
% clr = {'k','m'};
% nAn = size(norm_comp{2},2);

%% INPUTS
uni = unique({rawS.inf}); uni = fliplr(uni); % Unique infusions
nComp = length(uni); % Number of comparisons to be made
norm_comp = cell(1,nComp);
for x = 1:nComp
    norm_comp{x} = norm(:,[find(strcmp({rawS.inf}, uni{x}))]); end % Extract normalized FFT ouput for matching infusions into single matrix
lbl_inf = uni;
lbl_fp = rawS(1).fp_lbl;
switch lbl_fp; case 'ACh'; clr = {'k','g'}; case 'DA'; clr = {'k','m'}; end
nAn = size(norm_comp{1},2);

%% SAVE
%EXAMPLE: save('DataLocation\DataName.mat','norm_comp','f','flog','lbl_inf','lbl_fp');
% Must save variables: 'norm_comp','f','flog','lbl_inf','lbl_fp'

%% LOAD pre-analyzed signals into workspace
if ~exist('norm_comp','var')
    fPath = 'R:\tritsn01labspace\'; 
    [fName,fPath] = uigetfile([fPath,'*.mat'], 'Select data file that contains norm_comp f flog lbl_inf lbl_fp','MultiSelect','Off');
    if ~exist('norm_comp','var')
        error('No variable norm_comp in workspace.');
    end
end

%% LOAD fluorophore (GFP/tdTomato) signal
if ~exist('norm_gfp','var') && ~exist('norm_antd1d2','var') && ~exist('norm_tdt','var')
    [fName_flu,fPath_flu] = uigetfile([fPath,'*.mat'],'Select GFP+tdTomato for FFT file','MultiSelect','Off');
    load(fullfile(fPath_flu,fName_flu));
end
% fluorophore = menu('Fluorophore FFT signal to subtract','green','red');
switch lbl_fp
    case 'ACh'; sub = norm_gfp; % FFT ouput: GFP fluorescence signal, average over n = 3 mice
    % case 'DA'; sub_2 = norm_tdt; % FFT ouput: tdTomato fluorescence signal, average over n = 3 mice
    case 'DA'; sub = norm_antd1d2; % FFT ouput: rDA fluorescence signal during infusion of DA receptor antagonist, average over n = 4 mice
end
%% SUBTRACT stable fluorophore (GFP/tdTomato) signal
sub_comp = cell(1,nComp);
for y = 1:nComp
    for x = 1:size(norm_comp{y},2)
        sub_comp{y}(:,x) = norm_comp{y}(:,x) - nanmean(sub,2); 
    end
end

%% AUC in specified frequency band
range_auc = [0.5 4]; % Range for calculation of area under the curve, in Hz
auc = nan(nAn, nComp);
f_sub = 10.^(flog); % Regenerate frequency vector from log(freq)
r_auc = [find(f_sub == range_auc(1)):find(f_sub == range_auc(2))]; % extract AUC from [0.5 4] Hz
for y = 1:nComp
    for x = 1:nAn
        auc(x,y) = trapz(sub_comp{y}(r_auc,x))/length(r_auc);
    end
end

%% PLOT
ds = 50; % Downsample when plotting so figure is smaller in size
fig = figure; fig.Position(3) = 1000;
subplot(1,2,1); hold on
    for y = 1:nComp
        plotme = sub_comp{y}((1:ds:end),:);
        shadederrbar(flog(1:ds:end), nanmean(plotme,2), SEM(plotme,2), clr{y}); % Plot FFT averaged across all recordings
    end
    plot([-2 2],[0 0],'-','Color',[0 0 0 0.3]); % Line at 0 power
    plot([flog(f == range_auc(1)) flog(f == range_auc(1))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 0.5 Hz
    plot([flog(f == range_auc(2)) flog(f == range_auc(2))],[-0.1 0.4],'-','Color',[0 0 0 0.3]); % Line at 4 Hz 
    xlabel('Frequency'); ylabel('Power (a.u.)');
    xlim([-1 flog(f == 50)]); xticks([-2:2]); xticklabels({'0.01','0.1','1','10','100'});
    title(sprintf('FFT (n = %d): %s', nAn, lbl_fp)); axis square
subplot(1,2,2); hold on
    plot(auc', '--.k', 'MarkerSize', 20); 
    errorbar([0.75 2.25],nanmean(auc), SEM(auc,1), '.', 'MarkerSize', 20, 'Color', clr{2});
    xlim([0.5 2.5]); xticks([1 2]); xticklabels(lbl_inf);
    ylabel('Power (a.u.)'); ylim([0 0.5]); yticks([0:0.1:0.5]);
    [~,p] = ttest(auc(:,1),auc(:,2));
    title(sprintf('AUC p = %1.2f',p)); axis square
movegui(gcf,'center');