%% startup file for generalized Golub-Kahan codes, 
%   including gen-LSQR and gen-HyBR
%
% These codes were used to generate the figures in the paper:
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Chung & Saibaba (2017)

directory = pwd;
path(directory, path)

% RestoreTools
path([directory, '/RestoreTools/Classes'], path)
path([directory, '/RestoreTools/TestData'], path)

% Regularization Tools
path([directory, '/REGU'], path)

% Hybrid Codes
path([directory, '/genHyBR'], path)

% Toeplitz
path([directory, '/toeplitz'], path)

close all
clear all
clc