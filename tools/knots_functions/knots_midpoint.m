function [x,w]=knots_midpoint(n,x_a,x_b,whichrho)


% [x,w]=KNOTS_MIDPOINT(n,x_a,x_b) 
%
% implements the univariate n-points midpoint quadrature rule, i.e.,
% divides [x_a,x_b] in n subintervals of length h = (x_b - x_a)/n and returns as x the vector of midpoints
% The weights are all equal to  1/n, i.e., weight function rho(x) = 1/(x_b-x_a)
%
% [x,w]=KNOTS_MIDPOINT(n,x_a,x_b,'prob') 
%
% is the same as [x,w] = KNOTS_MIDPOINT(n,x_a,x_b)  above 
%
% [x,w]=KNOTS_MIDPOINT(n,x_a,x_b,'nonprob') 
%
% sets the weights to (x_b-x_a)/n,  i.e. weight function rho(x)=1


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

h = (x_b - x_a)/n;
x = linspace(x_a + h/2,  x_b - h/2, n);

if nargin==3
    whichrho = 'prob';
end

switch whichrho
    
    case 'nonprob'
        w=(x_b-x_a)/n*ones(1,n);

    case 'prob'
        w=1/n*ones(1,n);

    otherwise
    error('SparseGKit:WrongInput','4th input not recognized')
    
end
