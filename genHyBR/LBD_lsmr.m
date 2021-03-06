function [U, B, V, alpha, v] = LBD_lsmr(A, U, B, V, P_le, P_ri, options, alpha, v)
%
%     [U, B, V] = LBD_lsmr(A, U, B, V, P_le, P_ri, options)
%
%  Perform one step of Lanczos bidiagonalization with or without
%  reorthogonalization, WITHOUT preconditioner here.
%
% Input:
%          A - matrix
%       U, V - accumulation of vectors
%          B - bidiagonal matrix
% P_le, P_ri - input ignored
%    options - structure from HyBR (see HyBRset)
%
% Output:
%       U, V - updated "orthogonal" matrix
%          B - updated bidiagonal matrix
%
%  Refs:
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%
%   J.Chung and J. Nagy 3/2007

% Determine if we need to do reorthogonalization or not.
reorth = strcmp(HyBRget(options,'Reorth'), {'on'});

m = size(U,1);
k = size(B,2)+1;

if reorth % Need reorthogonalization
  if k == 1
    v = A'*U(:,k);
    alpha = norm(v);
    v = v / alpha;
  end
  u = A*v - alpha*U(:,k);
  
  for j = 1:k
    u = u - (U(:,j)'*u)*U(:,j);
  end
  
  beta = norm(u);
  u = u / beta;
  U = [U, u];
  
  V = [V, v];
  B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
  
  % For next iteration
  k = size(B,2)+1;
  v = A'*U(:,k) - B(k, k-1)*V(:,k-1);
  
  
  % Reorthogonalize
  for j = 1:k-1
    v = v - (V(:,j)'*v)*V(:,j);
  end
  alpha = norm(v);
  v = v / alpha;
  
  
  
  
  
else % Do not need reorthogonalization, save on storage
  if k == 1
    v = A'*U(:);
    alpha = norm(v);
    v = v / alpha;
  end
  
  u = A*v - alpha*U(:);
  
  beta = norm(u);
  u = u / beta;
  U = u(:);
  
  
  V = [V, v];
  B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
  
  % For next iteration
  k = size(B,2)+1;
  v = A'*U(:) - B(k, k-1)*V(:,k-1);
  alpha = norm(v);
  v = v / alpha;
end




