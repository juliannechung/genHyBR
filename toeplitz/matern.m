function k = matern(r,nu,ell)
%
%  Implementation of the Matern covariance kernel
%
%   Input: r is the radial distance
%         nu is a parameter that controls smoothness of the stochastic process
%         ell controls the correlation length
%   Output: k - Matern kernel
%
%   "Generalized Hybrid Iterative Methods for
%       Large-Scale Bayesian Inverse Problems"
%       - Chung and Saibaba, SISC, 2017
%
% Chung & Saibaba (2017)

scale = sqrt(2*nu)*r/ell;
fact  = 2.^(1-nu)./gamma(nu);

k = fact*(scale.^nu).*besselk(nu,scale);

%Fix the scaling issue at r = 0; K_nu(0) numerically evalates to inf
k(find(isnan(k))) = 1;

%For nu = inf, the square exponential covariance kernel
if nu == inf
  k = exp(-((r/ell).^2)/2);
end

end