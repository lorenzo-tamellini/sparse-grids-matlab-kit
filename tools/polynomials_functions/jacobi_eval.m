function L = jacobi_eval(x,k,alpha,beta,a,b)

% L = jacobi_eval(x,k,alpha,beta,a,b)
%
% returns the values of the k-th Jacobi polynomial orthoNORMAL in [a,b] w.r.t to 
% rho(x)=Gamma(alpha+beta+2)/(Gamma(alpha+1)*Gamma(beta+1)*(b-a)^(alpha+beta+1))*(x-a)^alpha*(b-x)^beta
% in the points x (x can be a matrix as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, 
% L_1(x) = (2*x-a-b)/(a-b)-(beta-alpha)/(alpha+beta+2)

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------



% this function expresses L as a function of the standard Jacobi "probabilistic" polynomial (i.e. orthoGONAL w.r.t. 
% rho(x)=Gamma(alpha+beta+2)/(Gamma(alpha+1)*Gamma(beta+1)*2^(alpha+beta+1))*(1+x)^alpha*(1-x)^beta,
% which are recursively calculated through the function standard_jac_eval, coded below in this .m file

% first compute the transformation of x (referred to Beta(alpha,beta,x_a,x_b) ) to z, the standard Beta(alpha,Beta,-1,1)

z = ( 2*x - a - b ) / ( a - b ) ;

% calculate the standard Legendre polynomials in t
L = standard_jac_eval(z,k,alpha,beta);
% modify L to take into account normalizations. 
if k>1
    C = gamma(k+alpha)*gamma(k+beta)/(gamma(alpha+1)*gamma(beta+1))*gamma(alpha+beta+2)/((2*(k-1)+alpha+beta+1)*(k+alpha+beta)*factorial(k-1));
    L = L / sqrt ( C );
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function L = standard_jac_eval(x,k,alpha,beta)

% L = standard_jac_eval(x,k,alpha,beta) 
%
% returns the values of the k-th standard Jacobi "probabilistic" polynomial (i.e. orthoGONAL w.r.t. 
% rho(x)=Gamma(alpha+beta+2)/(Gamma(alpha+1)*Gamma(beta+1)*2^(alpha+beta+1))*(1+x)^alpha*(1-x)^beta
% in the points x (x can be a vector as well)
%
% N.B. the polynomials start from k=0: L_0(x) = 1, L_1(x) = x-(beta-alpha)/(alpha+beta+2)

% base steps

% read this as L(k-2)
L_2=ones(size(x));

% and this as L(k-1)
L_1=x-(beta-alpha)/(alpha+beta+2);

if k==0
      L=L_2;
      return
elseif k==1
      L=L_1;
      return
else
      % recursive step
      for ric=2:k
            c1 = (beta^2-alpha^2)./((2*(ric-1)+alpha+beta).*(2*(ric-1)+alpha+beta+2));
            c2 = 4*(ric-1).*(ric-1+alpha).*(ric-1+beta).*(ric-1+alpha+beta)./((2*(ric-1)+alpha+beta).^2.*(2*(ric-1)+alpha+beta+1).*(2*(ric-1)+alpha+beta-1));
            L = (x-c1).* L_1 - c2*L_2;
            L_2=L_1;
            L_1=L;
      end
      return
end

