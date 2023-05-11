function L = generalized_lagu_eval_multidim(X,k,alpha,beta)

% L = GENERALIZED_LAGU_EVAL_MULTIDIM(X,k,alpha,beta)
%
% evaluates the multidim. generalized Laguerre polynomial of order k (multi-index) orthonormal on [0,+inf)^N 
% with respect to rho=prod_i
% beta_i^(alpha_i+1)/Gamma(alpha_i+1)*x^alpha_i*exp(-beta_i*x), alpha_i>-1,beta_i>0
% on the list of points X (each column is a point in R^N).
% Alpha and beta can be scalar values.

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------


[N, nb_pts] = size(X); % N here is the number of dimensions

% take care of the fact that alpha,beta may be scalar values

if length(alpha)==1
    alpha=alpha*ones(1,N);
    beta=beta*ones(1,N);    
end

% L is a product of N polynomials, one for each dim. We store the evaluation
% of each of these polynomials as rows of a matrix. We do not compute the 
% polynomials for k=0 though, we know already they will be 1

nonzero_n = find(k~=0);
if isempty(nonzero_n)
    L = ones(1,nb_pts);
else
    M = zeros(length(nonzero_n),nb_pts);
    j=0;
    for n=nonzero_n
        j=j+1;
        M(j,:) = generalized_lagu_eval(X(n,:),k(n),alpha(n),beta(n));
    end
    L = prod(M,1);
end
