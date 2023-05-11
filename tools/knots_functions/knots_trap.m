function [x,w]=knots_trap(n,x_a,x_b,whichrho)


% [x,w]=KNOTS_TRAP(n,x_a,x_b) 
%
% implements the univariate n-points trapezoidal quadrature rule with weight function rho(x) = 1/(x_b-x_a), 
% i.e., x = linspace(x_a,x_b,n) and w = 1/(x_b - x_a) * [h/2 h h ... h h/2], with h = (x_b - x_a)/(n-1) 
%
% [x,w]=KNOTS_TRAP(n,x_a,x_b,'prob') 
%
% is the same as [x,w] = KNOTS_TRAP(n,x_a,x_b)  above 
%
% [x,w]=KNOTS_TRAP(n,x_a,x_b,'nonprob') 
%
% sets the weights to w = [h/2 h h ... h h/2],  i.e. weight function rho(x)=1


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% if n<3
%     error('SparseGKit:WrongInput','at least 3 points needed for trap rule')
% end
% 
% h = (x_b - x_a)/(n-1);
% x = linspace(x_a,x_b,n);
% 
% if nargin==3
%     whichrho = 'prob';
% end
% 
% switch whichrho
%     
%     case 'nonprob'
%         w= [h/2 h*ones(1,n-2) h/2];
% 
%     case 'prob'
%         w= [h/2 h*ones(1,n-2) h/2]/(x_b - x_a);
% 
%     otherwise
%         error('SparseGKit:WrongInput','4th input not recognized')
%     
% end


if n == 1 % reduces to midpoint rule
    h = (x_b - x_a);
    x = (x_a+x_b)/2;
    w= h;
    
else
    h = (x_b - x_a)/(n-1);
    x = linspace(x_a,x_b,n);
    w= [h/2 h*ones(1,n-2) h/2];
end


if nargin==3
    whichrho = 'prob';
end


switch whichrho
    
    case 'nonprob'

    case 'prob'
        w= w/(x_b - x_a);

    otherwise
        error('SparseGKit:WrongInput','4th input not recognized')
    
end
