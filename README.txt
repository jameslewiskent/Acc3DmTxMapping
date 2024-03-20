This code generates and plots figures for 'Accelerated 3D Multi-Channel B1+ Mapping at 7T for the Brain and Heart'

Version 2.0

Install:
1. Please download the an unzip the code.
2. Start MatLab and navigate to the folder'Acc3DmTxMapping-main'
3. Open 'PlotSimulationFigures.m'
4. Run each section of code which will run simulations and produce plots for Figures 2, S3 and S4.

Notes: Some data is available as 'SyntheticBodyImages.mat' in the 'Data' directory. Some undersampling masks are also present, but should also be generated if they do not already exist. If a new mask is requested but produces an error, try going to the 'Required Functions\ESPIRiTvdPoisMex' directory and use the command 'mex vdPoisMex.c' in MatLab's Command Window to build mex.


Required MATLAB Toolboxes:
1. MATLAB (>=R2020b)
2. Signal Processing Toolbox
3. Image Processing Toolbox
4. Statistics and Machine Learning Toolbox
5. Parallel Computing Toolbox
6. MATLAB Parallel Server
7. Polyspace Bug Finder


Third party software downloaded from the internet was used extenstively in this project, these include:

SPIRiT V0.3 Toolbox - https://people.eecs.berkeley.edu/~mlustig/Software.html
1. M. Lustig and J. Pauly SPIRiT: Iterative Self-Consistent Parallel Imaging Reconstruction from Arbitrary k-Space Sampling MRM 2010;64(2):457-71
2. Zhang T, Pauly JM, Vasanawala SS, Lustig M. Coil Compression for Accelerated Imaging with Cartesian Sampling MRM 2013;69(2):571-82]
3. M. Uecker, P. Lai, M. J. Murphy, P. Virtue, M. Elad, J. M. Pauly, S. S. Vasanawala and M. Lustig, ESPIRiT – An Eigenvalue Approach to Autocalibrating Parallel MRI: Where SENSE meets GRAPPA, MRM 2013 published on-line
4. P. Shin, P.E.Z. Larson, M. A. Ohliger, M.Elad, J. M. Pauly, D. B. Vigneron and M. Lustig, Calibrationless Parallel Imaging Reconstruction Based on Structured Low-Rank Matrix Completion, 2013, accepted to Magn Reson Med.

TxLR - https://github.com/mchiew/txlr_paper
1. Hess AT, Dragonu I, Chiew M. Accelerated calibrationless parallel transmit mapping using joint transmit and receive low‐rank tensor completion. Magn Reson Med 2021;86:2454–2467 doi: 10.1002/mrm.28880.