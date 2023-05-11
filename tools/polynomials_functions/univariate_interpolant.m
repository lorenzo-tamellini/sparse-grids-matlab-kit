function f_values = univariate_interpolant(grid_points,function_on_grid,non_grid_points) 

% UNIVARIATE_INTERPOLANT interpolates a univariate function on grid points, i.e. evaluates the   
% Lagrange polynomial approximation on a generic point of the parameter space.
%   
% F_VALUES = UNIVARIATE_INTERPOLANT(GRID_POINTS,FUNCTION_ON_GRID,NON_GRID_POINTS) evaluates the
%       lagrangian interpolant of a vector function F: R -> R^V based on the points contained in the vector GRID_POINTS.
%       GRID POINTS is a row vector with the nodes of interpolation
%       FUNCTION_ON_GRID is a matrix containing the evaluation of F on the points GRID_POINTS. 
%       Its dimensions are V x length(GRID_POINTS)
%       NON_GRID_POINTS is a row vector of points where one wants to evaluate the polynomial approximation.
%       F_VALUES is a matrix containing the evaluation of the function F in each of the NON_GRID_POINTS. 
%       Its dimensions are V X length(NON_GRID_POINTS)


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------

[V,nb_grid_points] = size(function_on_grid);

if length(grid_points) ~= nb_grid_points
    error('SparseGKit:WrongInput','the number of function evaluations does not match the number of interpolation knots')
end

nb_pts = length(non_grid_points);

f_values = zeros(V,nb_pts);

nnn = 1:nb_grid_points;
    
for k=1:nb_grid_points
    f_values =  f_values + function_on_grid(:,k)*lagr_eval(grid_points(k), grid_points(nnn~=k),non_grid_points);
end 

