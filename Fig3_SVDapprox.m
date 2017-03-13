%% Fig3_SVDapprox.m
%
% This script investigates the approximation of the singular values of B_k
% and shows the interlacing property for \bar{B}_k
%
% This code was used to generate Figure 3 in the paper:
%
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Chung & Saibaba (2017)

%% Problem Setup
%   Generate Test Problem
n = 256;
[A, b, x_true] = heat(256); % From Per Christian Hansen's REGU toolbox

% Add noise to the data and solve the system
nlevel = .02;
R = speye(size(A,1));
[N,sigma] = WhiteNoise(b(:),nlevel);
bn = b(:) + N(:);

%% Select Q
covmat = 3; % Gausian 1D

%%%%%%%%%%%%%%%%%%%%%%1d test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = 0;
xmax = 1;

k = @(r) exp(-r./5.0);
theta = 1.0;

%Build row
Qr = createrow(xmin,xmax,n,k,theta);

%By the brute force approach
x = linspace(xmin,xmax,n);
[X,Y] = meshgrid(x,x);
Rvec = abs(X-Y)./theta;
Q = k(Rvec);

%% Get the singular values of Ahat
LQinv = chol(Q,'lower');
R_inv = diag(1./diag(R));
LR = diag(sqrt(1./diag(R)));
LRinv = diag(sqrt(diag(R)));
Ahat = LR*A*LQinv;

s_A = svd(full(A));
s_Ahat = svd(full(Ahat));

%% Now get the svd of B_k for various k
maxit = n;
input = HyBR_lsmrset('InSolv', 'none', 'x_true', x_true(:),'Iter', maxit,'Reorth','On');
[x, output] = genHyBR(A, bn(:), Q, R, input);
B = output.B;
V = output.V; U = output.U;

%% Include singular values of LSMR submatrix
figure, semilogy(s_Ahat,'r-'), hold on
for k = [5, 20, 35, 50]
  s = svd(B(1:k+1,1:k));
  plot(s,'k+-')
  plot(s.^2,'m*-')
  
  betabar = B(k+1,k)*B(k+1,k+1);
  Bsub = [B(1:k+1,1:k)'*B(1:k+1,1:k); zeros(1,k-1),betabar];
  s = svd(Bsub);
  plot(s,'bo-')
end
xlabel('k ')
axis([1,50,10^-12,50])

%% Include singular values of LSMR submatrix (NO REORTHOGONALIZATION)
maxit = n;
input = HyBR_lsmrset('InSolv', 'none', 'x_true', x_true(:),'Iter', maxit);
[x, output] = genHyBR(A, bn(:), Q, R, input);
B_no = output.B;
%%
figure,
semilogy(s_Ahat,'r-'), hold on
for k = [5, 20, 35, 50]
  s = svd(B_no(1:k+1,1:k));
  plot(s,'k+-')
  plot(s.^2,'m*-')
  
  betabar = B_no(k+1,k)*B_no(k+1,k+1);
  Bsub = [B_no(1:k+1,1:k)'*B_no(1:k+1,1:k); zeros(1,k-1),betabar];
  s = svd(Bsub);
  plot(s,'bo-')
end
h = legend('$\sigma_i(\hat{A})$', '$\sigma_i(B_k)$','$\sigma_i^2(B_k)$','$\sigma_i(\bar{B}_k)$','Location','sw');
set(h,'Interpreter','latex','FontSize',16)
xlabel('k')
axis([1,50,10^-12,50])

return
