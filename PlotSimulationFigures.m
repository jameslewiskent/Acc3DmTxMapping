%% Generate Data for and Plot Simulation Figures

% Version 1.0

% NOTES: 
% 1) Some of this script could take a while to run. If the masks are not
% available for the settings requested, the script will try to generate new 
% masks.  I recommend lowering N_Repeats to 5 to speed up code. 

% 2) If genCSmasks not running, go to 'Required
% Functions\ESPIRiTvdPoisMex' directory and use the command 'mex
% vdPoisMex.c' in MATLAB's Command Window to build mex
close all; clear all;
addpath(genpath('Required Functions')); addpath(genpath('Data'));

% Shim settings (do not change as these were used for generating synthetic images)
settings.Shim_Setting1 = exp((2*pi*1i*(1:8))/8); % CP+
settings.Shim_Setting2 = exp((2*pi*1i*2*(1:8))/8); % CP2+
%% First lets define some common variables (you can freely change any of these)
settings.sx = 24; settings.sy = 24; % Size to crop kspace down to
settings.fovx = 25; settings.fovy = 25;
settings.ellipse = 0; settings.pp = 0;
settings.accelerations = [1,2:2:10]; % Requested Acceleration Factors

settings.Mask_N = 1000; % Number of unique masks (1000 used in paper)
settings.NRepeats = 5; % Repeat with random channel-wise complex Gaussian noise and different undersampling masks (50 repeats used in paper, but requires minimum 500 masks be generated)

% Recon parameters
settings.sz = [32 32]; % Size to reconstruct to (Original Data is [139 178])
settings.niters_array = 50; % TxLR iterations
settings.pSNR = 60; % peak SNR
settings.Recon_Type = 'TxLR'; % Options include 'SENSE' 'GRAPPA' 'TxLR' 'sTxLR' or 'jTxLR' 
%% Generate data and plot Figure 2 (Also plots Supporting Figure 1)
settings.calib = 4; % calibration region
settings.type = 'square'; % aspect of figure

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(settings)

% Plot Figure 2
plotfigure2(settings)

%% Generate data and plot Supporting Figure 3
settings.calib = 0; % calibration region
settings.type = 'long'; % aspect of figure

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(settings)

% Plot Supporting Figure 3 (no calibration region; Compare to Figure 2)
plotfigure2(settings)

%% Generate data and plot Supporting Figure 4
settings.calib = 4; % calibration region
settings.niters_array = 2:2:50; % TxLR iterations (2:1:50 in paper)
settings.NRepeats = 5; % Repeat with random channel-wise complex Gaussian noise and different undersampling masks (50 repeats used in paper)

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(settings)

% Plot Supporting Figure 4
plotsupportingfigureS4(settings)

