%% Generate Data for and Plot Simulation Figures
% Version 2.0

% NOTES: 
% 1) Some of this script could take a while to run. If the masks are not
% available for the settings requested, the script will try to generate new 
% masks.  I recommend lowering N_Repeats to 5 to speed up code. 

% 2) If genCSmasks not running, go to 'Required
% Functions\ESPIRiTvdPoisMex' directory and use the command 'mex
% vdPoisMex.c' in MATLAB's Command Window to build mex
close all; clearvars;
addpath(genpath('Required Functions')); addpath(genpath('Data'));

settings.Enc_Scheme = 'B1TIAMO';
% If synthetic body images don't exist, generate them
if ~(exist(['Data',filesep,'SyntheticBodyImages.mat'],'file') == 2)
    GenerateSyntheticImages(settings.Enc_Scheme);
end
%% First lets define some common variables (you can freely change any of these)
settings.sx = 24; settings.sy = 24; % Size to crop kspace down to
settings.fovx = 25; settings.fovy = 25;
settings.ellipse = 0; settings.pp = 0;
settings.accelerations = [1,2:2:10]; % Requested Acceleration Factors

settings.Mask_N = 1000; % Number of unique masks (1000 used in paper)
settings.NRepeats = 50; % Repeat with random channel-wise complex Gaussian noise and different undersampling masks (50 repeats used in paper, but requires minimum 500 masks be generated)

% Recon parameters
settings.sz = [139 178]; % Size to reconstruct to (Original Data is [139 178])
settings.niters_array = 50; % TxLR iterations
settings.pSNR = 60; % peak SNR
%% Generate data and plot Figure 2 and Supporting Figure 1
settings.calib = 4; % calibration region
settings.type = 'square'; % aspect of figure

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(settings)

% Plot Figure 2
plotfigure2(settings)

%% Plot Figure S2
CompareCroppedMapsToGT;
%% Generate data and plot Supporting Figure 4
settings.calib = 4; % calibration region
settings.type = 'long'; % aspect of figure

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(settings)

% Plot Supporting Figure 4 (no calibration region; Compare to Figure 2)
plotfigure2(settings)

%% Generate data and plot Supporting Figure 4
settings.calib = 4; % calibration region
settings.niters_array = 2:2:50; % TxLR iterations (2:1:50 in paper)

% Generate Data
SimulateUndersamplingofSyntheticBodyImages(settings)

% Plot Supporting Figure 4
plotsupportingfigureS5(settings)

