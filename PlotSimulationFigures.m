%% Generate Data for and Plot Simulation Figures
% NOTES: 
% 1) You need to have downloaded the 'SyntheticBodyImages.mat' file and
% placed this within the 'Data' directory.
% 3) Some of this script could take a while to run. I recommend lowering
% Mask_N to 50 and N_Repeats to 10 if you have limited computing power.
% 2) If genCSmasks not running, go to 'Required Functions\ESPIRiTvdPoisMex' directory and
% use the command 'mex vdPoisMex.c' in MATLAB's Command Window to build mex
close all; clear all;
addpath(genpath('Required Functions')); addpath(genpath('Data'));

% Shim settings (do not change as these were used for generating synthetic images)
Shim_Setting1 = exp((2*pi*1i*(1:8))/8); % CP+
Shim_Setting2 = exp((2*pi*1i*2*(1:8))/8); % CP2+
%% First lets define some common variables (you can freely change any of these)
sx = 24; sy = 24; % Size to crop kspace down to
fovx = 25; fovy = 25;
ellipse = 0; pp = 0;
accelerations = [1,2:2:10]; % Requested Acceleration Factors

Mask_N = 1000; % Number of unique masks (1000 used in paper)
NRepeats = 50; % Repeat with random channel-wise complex Gaussian noise and different undersampling masks (50 repeats used in paper, but requires minimum 500 masks be generated)

% Recon parameters
sz = [32 32]; % Size to reconstruct to (Original Data is [139 178])
niters = 50; % TxLR iterations
pSNR = 60; % peak SNR
Recon_Type = 'TxLR'; % Options include 'SENSE' 'GRAPPA' 'TxLR' or 'jTxLR'
%% Generate data and plot Figure 2 (Also plots Supporting Figure 1)
calib = 4; % calibration region
[masks] = genCSmasks(sx, sy, fovx, fovy, calib, ellipse, pp, accelerations, Mask_N); % generate masks if neccesary

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(Recon_Type,sx,sy,sz,calib,niters,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

% Plot Figure 2
filename = ['Simulated_',Recon_Type,'Recon_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats),'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat'];
load(['Data',filesep,'Synthetic Body Simulation Results',filesep,'ReconData',filesep,filename],'Maps','Maps_acc');

plotfigure2(Maps,Maps_acc,sz,accelerations,Shim_Setting1,'square')

%% Generate data and plot Supporting Figure 3
calib = 0; % calibration region
[masks] = genCSmasks(sx, sy, fovx, fovy, calib, ellipse, pp, accelerations, Mask_N); % generate masks if neccesary

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(Recon_Type,sx,sy,sz,calib,niters,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

% Plot Supporting Figure 3 (no calibration region)
% (Used for comparison to Figure 2)
filename = ['Simulated_',Recon_Type,'Recon_sx',num2str(sx),'_sy',num2str(sy),'_calib',num2str(calib),'_niters',num2str(niters),'_Repeats',num2str(NRepeats),'_ReconSize',[num2str(sz(1)),num2str(sz(2))],'.mat'];
load(['Data',filesep,'Synthetic Body Simulation Results',filesep,'ReconData',filesep,filename],'Maps','Maps_acc');
plotfigure2(Maps,Maps_acc,sz,accelerations,Shim_Setting1,'square')

%% Generate data and plot Supporting Figure 4
calib = 4; % calibration region
[masks] = genCSmasks(sx, sy, fovx, fovy, calib, ellipse, pp, accelerations, Mask_N); % generate masks if neccesary
niters_array = 2:2:50; % TxLR iterations

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(Recon_Type,sx,sy,sz,calib,niters_array,NRepeats,pSNR,masks,accelerations,Shim_Setting1,Shim_Setting2)

% Plot Supporting Figure 4



