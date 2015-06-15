function H = herm_eval_multidim(X,k,mu,sig)

% H = HERM_EVAL_MULTIDIM(X,k)
%
% evaluates the multidim Hermite polynomial order k (multi-index) orthonormal on [-inf,+inf]^N 
% on the list of points X (each column is a point in R^N)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


[N, nb_pts] = size(X); % N here is the number of dimensions

H = zeros(1,nb_pts);

% H is a product of N polynomials, one for each dim. I store the evaluation
% of each of these polynomials as columns of a matrix

M = 0*X;

% take care of the fact that mu,sigma may be scalar values

if length(mu)==1
    mu=mu*ones(1,N);
    sig=sig*ones(1,N);
end



for n=1:N
    M(n,:) = herm_eval(X(n,:),k(n),mu(n),sig(n));
end

H = prod(M,1);