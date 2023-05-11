function [x,w]=knots_beta(n,alpha,beta,x_a,x_b)

% [x,w]=KNOTS_BETA(n,alpha,beta,x_a,x_b) 
% 
% calculates the collocation points (x) 
% and the weights (w) for the gaussian integration 
% w.r.t to the weight function 
% rho(x)=Gamma(alpha+beta+2)/(Gamma(alpha+1)*Gamma(beta+1)*(x_b-x_a)^(alpha+beta+1))*(x-x_a)^alpha*(x_b-x)^beta
% i.e. the density of a Beta random variable with range [x_a,x_b], alpha,beta>-1.

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

if (n==1) 
      % standard node
      x =(beta-alpha)/(alpha+beta+2);
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

% modifies points according to x_a and x_b (the weigths are unaffected)
x = ((x_a+x_b) - (x_b-x_a)*x) / 2;


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
