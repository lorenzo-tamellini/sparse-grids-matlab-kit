function [x,w]=knots_exponential(n,lambda)

% [x,w]=KNOTS_EXPONENTIAL(n,lambda) 
% 
% returns the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the weight function 
% rho(x)=lambda*exp( -lambda*x ) 
% i.e. the density of an exponential random variable 
% with rate parameter lambda

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------


if (n==1) 
      % the point (rescaled if needed) 
      x=1/lambda; 
      % the weight is 1:
      w=1;
      return
end

% calculates the values of the recursive relation
[a,b]=coeflagu(n); 

% builds the matrix
JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);

% calculates points and weights from eigenvalues / eigenvectors of JacM
[W,X]=eig(JacM); 
x=diag(X)'; 
w=W(1,:).^2;
[x,ind]=sort(x); %#ok<TRSRT>
w=w(ind);

% modifies points according to lambda (the weigths are unaffected)
x=x/lambda;



%----------------------------------------------------------------------
function [a, b] = coeflagu(n)

if (n <= 1)  
    disp(' n must be > 1 '); 
    return; 
end
a=zeros(n,1); 
b=zeros(n,1); 

a(1)=1; 
b(1)=1;
for k=2:n 
    a(k)=2*(k-1)+1; 
    b(k)=(k-1)^2; 
end
