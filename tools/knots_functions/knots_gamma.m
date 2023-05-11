function [x,w]=knots_gamma(n,alpha,beta)

% [x,w]=KNOTS_GAMMA(n,alpha,beta) 
% 
% returns the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the weight function 
% rho(x)= beta^(alpha+1)/Gamma(alpha+1)*x^alpha*exp(-beta*x)
% i.e. the density of a Gamma random variable 
% with range [0,inf), alpha>-1, beta>0

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------


if (n==1) 
      % the point  
      x=(alpha+1)/beta;
      % the weight is 1:
      w=1;
      return
end

% calculates the values of the recursive relation
[a,b]=coeflagu_generalized(n,alpha); 

% builds the matrix
JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);

% calculates points and weights from eigenvalues / eigenvectors of JacM
[W,X]=eig(JacM); 
x=diag(X)'; 
w=W(1,:).^2;
[x,ind]=sort(x); %#ok<TRSRT>
w=w(ind);

% rescales points 
x=x/beta;



%----------------------------------------------------------------------
function [a, b] = coeflagu_generalized(n,alpha)

if (n <= 1)
    disp(' n must be > 1 '); 
    return; 
end
a=zeros(n,1); 
b=zeros(n,1); 

a(1)=alpha+1; 
b(1)=gamma(1+alpha);

k=2:n;
a(k)=2*(k-1)+alpha+1; 
b(k)=(k-1).*(k-1+alpha); 


