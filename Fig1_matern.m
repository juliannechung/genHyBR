%% Fig1_matern.m
%
% This script investigates various Matern kernels and realizations.
%
% This code was used to generate Figure 1 in the paper:
%   "Generalized Hybrid Iterative Methods for 
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, 2017
%
% Chung & Saibaba (2017)

clear all
clc
close all

r = linspace(0,3,100);

%Typical Matern covariance kernels 
k12 = matern(r,0.5,1);
k32 = matern(r,1.5,1);
k52 = matern(r,2.5,1);

figure,
plot(r,k12,'-',r,k32,'--',r,k52,':','LineWidth',2.0);
l = legend('$\nu = 1/2$','$\nu = 3/2$','$\nu = 5/2$');
set(l,'Interpreter','LaTex','fontsize',12)

%% Plot realizations
x = linspace(0,3,100); 
[X,Y] = meshgrid(x,x);
R = abs(X-Y);
Q12 = matern(R,0.5,1);
L = chol(Q12,'lower');
real12 = L*randn(100,1);
Q32 = matern(R,1.5,1);
L = chol(Q32,'lower');
real32 = L*randn(100,1);
Q52 = matern(R,2.5,1);
L = chol(Q52,'lower');
real52 = L*randn(100,1);

figure,
plot(r,real12,'-',r,real32,'--',r,real52,':','LineWidth',2.0);


