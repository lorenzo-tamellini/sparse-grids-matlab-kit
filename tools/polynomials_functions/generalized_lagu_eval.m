function L = generalized_lagu_eval(x,k,alpha,beta)

% L = generalized_lagu_eval(x,k,alpha,beta)
%
% returns the values of the k-th generalized Laguerre polynomial orthoNORMAL in [0,+inf) w.r.t to 
% rho=beta^(alpha+1)/Gamma(alpha+1)*x^alpha*exp(-beta*x) in the points x (x can be a matrix as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = - x*beta + alpha + 1

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------



% this function expresses L as a function of the standard generalized
% Laguerre "probabilistic" polynomial (i.e. orthoGONAL w.r.t. rho=x^alpha*exp(-x)/Gamma(alpha+1),
% which are recursively calculated through the function standard_generalized_lagu_eval, coded below in this .m file

% first compute the transformation of x (referred to Gamma(alpha,beta)) to z, the standard Gamma(alpha,1)

z = x * beta ;

% calculate the standard generalized Laguerre polynomials in z
L = standard_generalized_lagu_eval(z,k,alpha);
% modify L to take into account normalizations. 
if k>=1
    L = L / sqrt ( gamma(k+alpha+1)/(gamma(alpha+1)*factorial(k)) );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function L = standard_generalized_lagu_eval(x,k,alpha)

% L = standard_generalized_lagu_eval(x,k,alpha) 
%
% returns the values of the k-th standard generalized Laguerre "probabilistic" polynomial (i.e. orthoGONAL w.r.t. 
% rho=x^alpha*exp(-x)/Gamma(alpha+1) in the points x (x can be a vector as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = -x+alpha+1

% base steps

% read this as L(k-2)
L_2=ones(size(x));

% and this as L(k-1)
L_1= - x + alpha + 1;

if k==0
      L=L_2;
      return
elseif k==1
      L=L_1;
      return
else
      % recursive step
      for ric=2:k
            L = (- x + 2*(ric-1) + alpha + 1)/ric .* L_1 - (ric-1+alpha)/ric *L_2;
            L_2=L_1;
            L_1=L;
      end
      return
end


