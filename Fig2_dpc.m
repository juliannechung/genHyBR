%% Fig2_dpc.m
%
% This script investigates the Discrete Picard condition corresponding to
% the generalized SVD.  
%
% This code was used to generate Figure 2 in the paper:
%
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Chung & Saibaba (2017)

close all

%% Problem Setup
n = 256;
[A, b, x_true] = heat(n); % From Per Christian Hansen's REGU toolbox

% Add noise to the data
noisetype = 0; % White Noise
noise  = randn(size(b(:)));
sigma = 0.0001*(norm(b)/norm(noise));
N = sigma*noise;
b = b(:) + N(:);

R = speye(size(A,1));    
%% Discrete Picard plot for standard SVD
[U_0,S,V_0] = svd(A);
s = diag(S);
bhat_0 = U_0'*b;

%% Select Q
% covmat = 0; % Identitiy
covmat = 1; % Gausian
switch covmat
  case 0 % Identity
    Q = eye(n);
  case 1 %
    %%%%%%%%%%%%%%%%%%%%%%1d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xmin = 0;
    xmax = 1;
 
    k = @(r) exp(-10*r./5.0);
    theta = 1.0;

    %Build row
    Qr = createrow(xmin,xmax,n,k,theta);

    %By the brute force approach
    x = linspace(xmin,xmax,n);
    [X,Y] = meshgrid(x,x);
    Rvec = abs(X-Y)./theta;
    Q = k(Rvec);
end

%% Get the singular values of Ahat
LQinv = chol(Q,'lower');
R_inv = diag(1./diag(R));
LR = diag(sqrt(1./diag(R)));
LRinv = diag(sqrt(diag(R)));
Ahat = LR*A*LQinv;

[U,s_Ahat,V] = svd(full(Ahat));
s_Ahat = diag(s_Ahat);
Ubar = LRinv*U;
Vbar = LQinv*V;
bhat = Ubar'*(R_inv*b);

%% Put both Picard plots together
figure, semilogy(s,'b+-'), hold on
semilogy(s_Ahat,'m--')
semilogy(abs(bhat_0),'r.-', 'MarkerSize',8),
semilogy(abs(bhat),':','color',[0.5,0.5,0.5],'LineWidth',2.0),
semilogy(abs(bhat_0./s),'cs-'),
semilogy(abs(bhat)./s_Ahat,'ko-'),
axis tight

xlabel('j')
h = legend('$\sigma_j(A)$', '$\hat{\sigma_j}(\hat{A})$',...
  '$|u_j^T b|$', '$|\bar{u}_j^T b|$', ...
  '$|u_j^T b|/\sigma_j$','$|\bar{u}_j^T b|/\hat{\sigma_j}$');
set(h, 'Interpreter', 'latex','FontSize', 18,'Location','northeastoutside')
axis([0,200,10^-9,10^2])

return
