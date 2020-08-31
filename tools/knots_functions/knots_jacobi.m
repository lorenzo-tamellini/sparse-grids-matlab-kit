function [x,w]=knots_jacobi(n,alpha,beta,whichrho)

% [x,w]=KNOTS_JACOBI(n,alpha,beta) 
% 
% calculates the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the weight function 
% rho(x)=Gamma(alpha+beta+2)/(Gamma(alpha+1)*Gamma(beta+1)) x^alpha*(1-x)^beta
% i.e. the density of a Beta random variable with range [0,1], alpha,beta>-1.
%
%
% [x,w]=KNOTS_JACOBI(n,alpha,beta,'prob')  
% 
% is the same as [x,w]=KNOTS_JACOBI(n,alpha,beta) above
% 
%
% [x,w]=KNOTS_JACOBI(n,alpha,beta,'nonprob') 
%
% calculates the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the Jacobi weight function 
% rho(x)= (1-x)^alpha(1+x)^beta 
 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if nargin==3
    whichrho='prob';
end

if (n==1) 
      % standard node: t = (beta-alpha)/(alpha+beta+2)
      x=(beta-alpha)/(alpha+beta+2);
      % the weight is 1:
      w=1;
else
    
    % calculates the values of the recursive relation
    [a,b]=coefjac(n,alpha,beta); 

    % builds the matrix
    JacM=diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);

    % calculates points and weights from eigenvalues / eigenvectors of JacM
    [W,X]=eig(JacM); 
    x=diag(X)'; 
    w=W(1,:).^2;
    [x,ind]=sort(x); %#ok<TRSRT>
    w=w(ind);

end

if strcmp(whichrho,'prob') 
    
    % modifies points (the weigths are unaffected)
    x=(1-x)/2;
    
end


%----------------------------------------------------------------------
function [a, b] = coefjac(n,alpha,beta)

if (n <= 1),  
    disp(' n must be > 1 '); 
    return; 
end
a=zeros(n,1); 
b=zeros(n,1); 

a(1)=(beta-alpha)/(alpha+beta+2);
b(1)=2^(alpha+beta+1)*gamma(alpha+1)*gamma(beta+1)/gamma(alpha+beta+2);
a(2)=(beta^2-alpha^2)./((2+alpha+beta).*(2+alpha+beta+2)); 
b(2)=4*(alpha+1)*(beta+1)/((2+alpha+beta)^2*(alpha+beta+3));

k=3:n;
a(k)=(beta^2-alpha^2)./((2*(k-1)+alpha+beta).*(2*(k-1)+alpha+beta+2));
b(k)=(4*(k-1).*(k-1+alpha).*(k-1+beta).*(k-1+alpha+beta))./((2*(k-1)+alpha+beta).^2.*(2*(k-1)+alpha+beta+1).*(2*(k-1)+alpha+beta-1)); 
