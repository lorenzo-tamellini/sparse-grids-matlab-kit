function Hessian = hessian_sparse_grid(S,Sr,values_on_grid,domain,eval_point,h)
  

% HESSIAN_SPARSE_GRID computes the hessian of a scalar-valued function f: R^N -> R 
% by finite differences formulas applied to the sparse grids approximation of f.
% Contrary to derive_sparse_grid, this function can be evaluated at *a single point*. 
% 
% HESSIAN = HESSIAN_SPARSE_GRID(S,SR,VALUES_ON_GRID,DOMAIN,EVAL_POINT) computes the hessian of f
%           by finite differences (see below for default value of increment size).
%
%          S is a sparse grid struct, SR is the reduced version of S, VALUES_ON_GRID are the values of the interpolated
%          function on SR. 
%          DOMAIN is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] defining the lower and upper bound of the 
%          hyper-rectangle on which the sparse grid is defined. The finite differences increment size is chosen according to
%          to the length of each interval [an bn] as h_n = (b_n - a_n)/1E5 
%
%          EVAL_POINT is the point where the derivative must be evaluated. It is a column vector point, following the convention of the package
%
%          The output HESSIAN contains is the Hessian matrix, H(i,j) = \partial_{y_i} \partial_{y_j} sparse grid 
%
%
% HESSIAN = DERIVE_SPARSE_GRID(S,SR,VALUES_ON_GRID,DOMAIN,EVAL_POINTS,H) uses the input H as finite differences 
%           increment. H can be a scalar or a vector, in which case the n-th entry will be used as increment to approximate
%           the n-th component of the gradient


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% get dimensions
N=size(domain,2);   % the sparse grid is defined over an N-variate hypercube


Hessian = zeros(N,N);

switch nargin
    case 5
        h = zeros(1,N);
        a = domain(1,:);
        b = domain(2,:);
        
        for n=1:N
            h(n) = ( b(n)-a(n) ) / 1E5;
        end
    case 6
        if length(h) ==1
            h = h*ones(1,N);
        end
end

% the evaluation of the sparse grids at the point (i.e. the center of the stencil) can be done once and for all
f_0 = interpolate_on_sparse_grid(S,Sr,values_on_grid,eval_point);

for j=1:N
    
    %-------------------------------------------------------
    % the extra-diagonal part of the matrix (upper part)    
    for k = j+1:N

        % the formula for D_jk uses 4 points: (only values of j and k index are specified)
        % ( f_{j+1,k+1} - f_{j-1,k+1} - f_{j+1,k-1} + f_{j-1,k-1} ) / (4 h_j h_k)
        
        epsi = zeros(N,4); 
        
        % the first point is beta + he_j + he_k
        epsi([j k],1) = [ h(j)  h(k)];
        % the second point is beta - he_j + he_k
        epsi([j k],2) = [-h(j)  h(k)];
        % the third point is beta + he_j - he_k
        epsi([j k],3) = [ h(j) -h(k)];
        % the fourth point is beta -he_j - he_k
        epsi([j k],4) = [-h(j) -h(k)];

        % the evaluation
        f_evals = interpolate_on_sparse_grid(S,Sr,values_on_grid,eval_point+epsi);
        % the approx of Hessian entry
        Hessian(j,k) = dot(f_evals,[1 -1 -1 1]) / ( 4*h(j)*h(k) ); 
        
    end

    %-------------------------------------------------------
    % the extra-diagonal part of the matrix (lower part) is copied from the upper part that we have already computed
    for k = 1:j-1
        Hessian(j,k) = Hessian(k,j);        
    end
    
    %-------------------------------------------------------
    % then the diagonal entry
    % the formula for D_jj uses 3 points: (only values of j and k index are specified)
    % ( f_{j+1} - 2*f_{j} + f_{j-1} ) / h_j^2
    %
    % note that we already computed f_j once and for all before the loop (variable f_0)
    
    epsi = zeros(N,2);
    
    % the first point is beta + he_j 
    epsi(j,1) = h(j);
    % the second point is beta - he_j
    epsi(j,2) = -h(j);
    
    % the evaluation
    f_evals = interpolate_on_sparse_grid(S,Sr,values_on_grid,eval_point+epsi);
    % the approx of Hessian entry
    Hessian(j,j) = ( sum(f_evals) - 2*f_0) / h(j)^2;

    
    
end


