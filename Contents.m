%% Contents.m
%   This folder contains MATLAB code to accompany the paper:
%
%   "Generalized Hybrid Iterative Methods for
%                   Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
%   The DEMO codes require the following packages:
%       REGU Regularization Tools by Per Christian Hansen
%             http://www.imm.dtu.dk/~pcha/Regutools/
%       RestoreTools by James Nagy
%             http://www.mathcs.emory.edu/~nagy/RestoreTools/
% 
%   First run startup_genGK.m for setting paths
%
% Chung & Saibaba (2017)

%% Script files used to generate some figures in the paper
%
%   Fig1_matern.m          Generates kernels and realizations from the
%                             Matern covariance class. 
%
%   Fig2_dpc.m             Generates discrete Picard plots for SVD and GSVD.
%
%   Fig3_SVDapprox.m       Generates plots illustrating the approximation of
%                             the singular values of B_k and shows the
%                             interlacing property for \bar{B}_k.
%

%% DEMOs on how to use genHyBR
%   DEMO_tomo.m           Sets up and runs a 2D tomgraphy problem
%
%   DEMO_deblur.m         Sets up and runs a 2D image deblurring problem
%


%% Supporting MATLAB functions
%
% genHyBR.m             Main code to run generalized HyBR
%                       Syntax: [x_out, output] = genHyBR(A, b, Q, R, options)

% genHyBR_lsmr.m        Main code to run generalized HyBR - LSMR version
%                       Syntax: [x_out, output] = genHyBR_lsmr(A, b, Q, R, options)
