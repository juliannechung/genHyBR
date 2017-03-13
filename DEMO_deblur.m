% DEMO_deblur.m
%
% This script sets up and computes reconstructions for a 2D deblurring problem.
%
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Note: This code requires RestoreTools: 
%             http://www.mathcs.emory.edu/~nagy/RestoreTools/
%
% Chung & Saibaba (2017)

close all

n = 256;
load Grain
A = psfMatrix(PSF,[129,129],'reflexive');
b = A*x_true;

% Add noise to the data
noise  = randn(size(b(:)));
sigma = 0.05*(norm(b)/norm(noise));
N = sigma*noise;
bn = b(:) + N(:);
figure, imshow(reshape(bn,n,n),[])

R = speye(size(A,1));

%% Select Q
%Example application of Qx in 2D
%%%%%%%%%%%%%%%%%%%%%%2d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider the box (0,1)^2
xmin = [0 0];         %Coordinates of left corner
xmax = [1 1];         %Coordinates of right corner
nvec = [n, n];        %Number of points in grid

nu = 1.5; ell = .007;
k = @(r) matern(r,nu,ell);
%Additional parameters governing length scales. Can be used to model
%nonisotropic covariance kernels
theta = [1.0 1.0];      %For now set them as isotropic

% Build row/column of the Toeplitz matrix
Qr = createrow(xmin,xmax,nvec,k,theta);
Qfun = @(x)toeplitzproduct(x, Qr, nvec);
Q = funMat(Qfun,Qfun, nvec.^2);

%% Compare different methods to get a reconstruction
maxit = 100;

% LSQR
input = HyBRset('InSolv', 'none', 'x_true', x_true(:),'Iter', maxit);
[x_cg_n, output_lsqr] = HyBR(A, bn(:), [],input);

% HyBR with optimal regularization parameter
solver = 'tikhonov';
input = HyBRset('InSolv', solver,'RegPar','optimal', 'x_true', x_true(:),'Iter', maxit);
[x_hybr, output_hybr] = HyBR(A, bn(:), [],input);
fprintf('HyBR alpha optimal: %.4e\n',output_hybr.alpha)

% gen-LSQR - no regularization
input = HyBR_lsmrset('InSolv', 'none','x_true', x_true(:),'Iter', maxit);
[x_genlsqr, output_genlsqr] = genHyBR(A, bn(:), Q, R, input);

% gen-HyBR with optimal regularization parameter
input = HyBR_lsmrset('InSolv', solver,'RegPar','optimal', 'x_true', x_true(:),'Iter', maxit);
[x_genHyBR_opt, output_genHyBR_opt] = genHyBR(A, bn(:), Q, R, input);
fprintf('genHyBR alpha optimal: %.4e\n',output_genHyBR_opt.alpha)


%% Plot relative errors
figure, 
plot([1; output_lsqr.Enrm], 'b--*', 'LineWidth', 1),hold on
plot([1; output_hybr.Enrm],'b-', 'LineWidth', 1)
plot(output_genlsqr.Enrm, 'r-x', 'LineWidth', 1),
plot( output_genHyBR_opt.Enrm, 'r-.', 'LineWidth', 1),

ylabel('relative error')
xlabel('iteration')
l = legend('LSQR','HyBR-opt', 'gen-LSQR', 'gen-HyBR-opt');
axis([1,maxit,.2,.6])

%% Image Reconstructions
figure,
subplot(2,2,1), imshow(reshape(x_true,n,n),[]), title('True')
subplot(2,2,2), imshow(reshape(x_cg_n,n,n),[]), title('LSQR')
subplot(2,2,3), imshow(reshape(x_hybr,n,n),[]), title('HyBR-opt')
subplot(2,2,4), imshow(reshape(x_genHyBR_opt,n,n),[]), title('genHyBR-opt')


