function L = jacobi_prob_eval_multidim(X,k,alpha,beta,a,b)

% L = JACOBI_PROB_EVAL_MULTIDIM(X,k,alpha,beta,a,b)
%
% evaluates the multidim. probabilistic Jacobi polynomial of order k
% (multi-index) orthonormal on [a1,b1] x [a2,b2] x [a3,b3] x ... [aN,bN] 
% with respect to 
% rho=prod_i Gamma(alpha_i+beta_i+2)/(Gamma(alpha_i+1)*Gamma(beta_i+1)*(b-a)^(alpha_i+beta_i+1))*(x-a)^alpha_i*(b-x)^beta_i,
% alpha_i,beta_i>-1 on the list of points X (each column is a point in R^N).
% Alpha, beta, a, and b can be scalar values.

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------


[N, nb_pts] = size(X); % N here is the number of dimensions

% take care of the fact that alpha,beta may be scalar values

if length(alpha)==1
    alpha = alpha*ones(1,N);
    beta = beta*ones(1,N);
end

if length(a)==1
    a = a*ones(1,N);
    b = b*ones(1,N);
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
        M(j,:) = jacobi_prob_eval(X(n,:),k(n),alpha(n),beta(n),a(n),b(n));
    end
    L = prod(M,1);
end
