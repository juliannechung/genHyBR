function [x, alpha] = Tikhonovsolver(U, S, V, b, options, B, beta, V_gk, m)
%
%         [x, alpha] = Tikhonovsolver(U, S, V, b, options, B, beta, V_gk, m)
%
%  This function computes a Tikhonov regularized LS solution to the
%  PROJECTED problem, using the identity matrix as the regularization operator.
%
%  Input:
%    U, S, V - SVD of A, where S is a column vector containing
%              singular (or spectral) values of A.
%          b - noisy image
%    options - structure with the following necessary fields:
%         RegPar: [value | GCV | {WGCV} | optimal] where
%                     (a) a value for the regularization parameter
%                     (b) a method for choosing a reg. parameter
%                                (i) 'GCV' - standard GCV
%                                (ii) 'WGCV' - weighted GCV (default)
%                     (c) 'optimal' - finds the optimal reg. parameter
%                                (requires x_true)
%         Omega: if RegPar is 'WGCV', then omega must be
%                           [value | {adapt}]
%                     (a) a value for omega
%                     (b) 'adapt' - use the adaptive method (default)
%            Vx: V_k * x_true (where V_k comes from Lanczos bidiagonlization)
%                 This vector is needed to compute the optimal parameter.
%          B - bidiagonal matrix
%       beta - norm of original rhs
%       V_gk - Golub-Kahan basis vectors (vectors multiplied with Q for gen-HyBR)
%          m - size of original problem
%
%  Output:
%           x - Tikhonov solution
%       alpha - regularization parameter
%
% J.Chung and J. Nagy 3/2007

bhat = U'*b;
bhat = bhat(:);

% Get values from the options structure.
alpha = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');

switch alpha
  case 'dp'
    nLevel = HyBRget(options,'nLevel',[],'fast');
    alpha = fminbnd('TikDP', 0, S(1), [], V, bhat, S, nLevel, B, beta, m);
    
  case 'gcv'
    alpha = fminbnd('TikGCV', 0, S(1), [], bhat, S);
    
  case 'wgcv'
    alpha = fminbnd('TikGCV', 0, S(1), [], bhat, S, omega);
    
  case 'upre'
    nLevel = HyBRget(options,'nLevel',[],'fast');
    alpha = fminbnd('TikUPRE', 0, S(1), [], bhat, S, nLevel);
    
  case 'optimal'
    x_true = HyBRget(options,'x_true',[],'fast');
    errhan = @(lambda)TikOPT(lambda, V, bhat, S, V_gk, x_true);
    options.TolX = eps;
    alpha = fminbnd(errhan, 0, S(1), options);
    
  otherwise
    if ~isa(alpha, 'double')
      error('Incorrect input argument for alpha.')
    end
end

D = abs(S).^2 + alpha^2;
bhat = conj(S) .* bhat(1:length(S));
xhat = bhat ./ D;
x = V * xhat;