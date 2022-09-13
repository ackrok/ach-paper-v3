%% Generating cannula output structure
%
% Anya Krok, September 2022
%
%% INPUTS
infOptions = {'saline','iGluR antag','nAChR antag','mAChR antag','D1/2R antag'};
infWindow = [20 40]; % window to analyze after infusion, in minutes
if ~exist('cannula'); cannula = struct; end

%% Select .mat files you want to add to summary data structu
fPath = 'R:\tritsn01labspace\'; 
[fName,fPath] = uigetfile([fPath,'*.mat'],'MultiSelect','On');
beh = extractBeh(fPath, fName); % extract data
inf = menu('Select infusion',infOptions); % select data type

cannula(inf).s = beh; % load into structure
cannula(inf).inf = infOptions(inf); % infusion
cannula(inf).win = infWindow.*ones(length(beh),2); % window to analyze after infusion, in minutes
