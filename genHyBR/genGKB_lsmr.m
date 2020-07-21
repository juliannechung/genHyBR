function [U, B, V, QV, alpha, v] = genGKB_lsmr(A, Q, R, U, B, V, QV, options, alpha, v)
%
%     [U, B, V, QV, alpha, v] = genGKB_lsmr(A, Q, R, U, B, V, QV, options, alpha, v)
%
%  Perform one step of generalized Golub-Kahan bidiagonalization for LSMR.
%
% Input:
%          A - matrix
%       Q, R - covariance matrices
%   U, V, QV - accumulation of vectors
%          B - bidiagonal matrix
%    options - structure from HyBR (see HyBRset)
%     alpha, v - additional pieces for LSMR
%
% Output:
%       U, V - updated "orthogonal" matrix
%          B - updated bidiagonal matrix
%   alpha, v - additional pieces for LSMR
%         QV - matrix of vectors Q*v_j
%
%  Refs:
%
%   J.Chung and A.Saibaba, 2017

% Determine if we need to do reorthogonalization or not.
reorth = strcmp(HyBR_lsmrget(options,'Reorth'), {'on'});

m = size(U,1);
k = size(B,2)+1;

if reorth % Need reorthogonalization
  if k == 1
    v = A'*(R\U(:,k));
    alpha = normM(v,Q);
    v = v / alpha;
  end
  
  Qv = Q*v;
  u = A*Qv - alpha*U(:,k);
  
  % Reorthogonalize U
  for j = 1:k
    temp = (U(:,j)'*(R\u));
    u = u - temp*U(:,j);
  end
  
  beta = normM(u, @(x)R\x);
  u = u / beta;
  U = [U, u];
  
  V = [V, v];
  B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
  
  % LSMR only: For next iteration
  k = size(B,2)+1;
  v = A'*(R\U(:)) - B(k, k-1)*V(:,k-1);
  
  % Reorthogonalize V
  for j = 1:k-1
    temp = (V(:,j)'*(Q*v));
    v = v - temp*V(:,j);
  end
  alpha = normM(v,Q);
  v = v / alpha;
 
else % Do not need reorthogonalization, save on storage
  if k == 1
    v = A'*(R\U(:));
    alpha = normM(v,Q);
    v = v / alpha;
  end
  Qv = Q*v;
  u = A*Qv - alpha*U(:);
  
  beta = normM(u, @(x)R\x);
  u = u / beta;
  U = u(:);
  
  V = [V, v];
  B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
  
  % LSMR only: For next iteration
  k = size(B,2)+1;
  v = A'*(R\U(:)) - B(k, k-1)*V(:,k-1);
  alpha = normM(v,Q);
  v = v / alpha;
  
end

QV = [QV, Qv];

end

function nrm = normM(v, M)
if isa(M, 'function_handle')
  Mv = M(v);
else
  Mv = M*v;
end
nrm = sqrt(v'*Mv);
end