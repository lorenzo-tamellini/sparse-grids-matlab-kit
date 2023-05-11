function L = lagu_eval(x,k,lambda)

% L = lagu_eval(x,k,lambda)
%
% returns the values of the k-th Laguerre polynomial orthoNORMAL in [0,+inf) w.r.t to rho=lambda*exp(-lambda*x), lambda>0, 
% in the points x (x can be a matrix as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = - x*lambda + 1 

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

% this function expresses L as a function of the standard Laguerre "probabilistic" polynomial (i.e. orthoGONAL w.r.t. rho=e^(-x),
% which are recursively calculated through the function standard_lagu_eval, coded below in this .m file

% first compute the transformation of x (referred to Exp(lambda)) to z, the standard Exp(1)

z = x * lambda ;

% calculate the standard Laguerre polynomials in z
L = standard_lagu_eval(z,k);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function L = standard_lagu_eval(x,k)

% L = standard_lagu_eval(x,k) 
%
% returns the values of the k-th standard Laguerre "probabilistic" polynomial (i.e. orthoGONAL w.r.t. rho=exp(-x), in the points x
% ( x can be a vector as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = -x+1

% base steps

% read this as L(k-2)
L_2=ones(size(x));

% and this as L(k-1)
L_1= -x+1;

if k==0
      L=L_2;
      return
elseif k==1
      L=L_1;
      return
else
      % recursive step
      for ric=2:k
            L = (- x + 2*(ric-1) + 1)/ric .* L_1 - (ric-1)/ric * L_2;
            L_2=L_1;
            L_1=L;
      end
      return
end


