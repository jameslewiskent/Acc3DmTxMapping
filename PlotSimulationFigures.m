%% Generate Data for and Plot Simulation Figures
% NOTE: You need to have downloaded the 'SyntheticBodyImages.mat' file and
% placed this within the 'Data' directory.
close all; clear all;
addpath(genpath('Acc3DmTxMapping'))

% Shim settings (do not change as these were used for generating synthetic images)
Shim_Setting1 = exp((2*pi*1i*(1:8))/8); % CP+
Shim_Setting2 = exp((2*pi*1i*2*(1:8))/8); % CP2+
%% First lets generate some masks
% Note: if not running, go to 'Required Functions\ESPIRiTvdPoisMex' directory and
% use the command 'mex vdPoisMex.c' in MATLAB's Command Window to build mex
sx = 24; sy = 24; % size to crop kspace down to
fovx = 25; fovy = 25;
calib = 0; % calibration region
ellipse = 0; pp = 0;
accelerations = [1,2:2:16]; % Requested Acceleration Factors
Mask_N = 100; % Number of unique masks (1000 used in paper, but takes a while so try less)
[masks] = genCSmasks(sx, sy, fovx, fovy, calib, ellipse, pp, accelerations, Mask_N);

%% Generate data for Figure 2 (Also plots Supporting Figure 1)
% These values can be changed as desired
sz = [32 32]; % Size to reconstruct to (Original Data is [139 178])
niters = 50; % TxLR iterations
NRepeats = 2; % Repeat with different complex channelwise Gaussian noise and different undersampling masks (50 repeats used in paper, but requires minimum 500 masks be generated)
pSNR = 60; % peak SNR

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(sx,sy,sz,calib,niters,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

%% Plot Figure 2
filename = ['Simulated_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats),'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat'];
load(['Data',filesep,'Synthetic Body Simulation Results',filesep,'ReconData',filesep,filename],'Maps','Maps_acc');
plotsyntheticheartmaps(Maps,Maps_acc,sz,accelerations,Shim_Setting1)

%% Generate data for Supporting Figure 3
calib = 0; % calibration region
load(['Data',filesep,'Undersampling Masks',filesep,'x',num2str(sx),'y',num2str(sy),'calib',num2str(calib),filesep,'masks.mat'],'masks')

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(sx,sy,sz,calib,niters,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

%% Plot Supporting Figure 3 (no calibration region)
% (Used for comparison to Figure 2)
filename = ['Simulated_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats),'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat'];
load(['Data',filesep,'Synthetic Body Simulation Results',filesep,'ReconData',filesep,filename],'Maps','Maps_acc');
plotsyntheticheartmaps(Maps,Maps_acc,sz,accelerations,Shim_Setting1)


%% Generate data for Supporting Figure 4
sx = 24; sy = 24; % Crop down to sx * sy
calib = 4; % calibration region
load(['Data',filesep,'Undersampling Masks',filesep,'x',num2str(sx),'y',num2str(sy),'calib',num2str(calib),'/masks.mat'],'masks')
niters_array = 2:2:50; % TxLR iterations

sz = [32 32]; % Size to reconstruct to (Original Data is [139 178])
niters = 50; % TxLR iterations
NRepeats = 10; % Repeat with different complex channelwise Gaussian noise and different undersampling masks
pSNR = 60; % peak SNR

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(sx,sy,sz,calib,niters_array,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

%% Plot Supporting Figure 4


