function res = quadrature_on_sparse_grid(f,S,map,weights_fact)

% res = QUADRATURE_ON_SPARSE_GRID(f,S)
%
% compute the integral of f using the sparse grid structure S.
%
% -> S can be either reduced or not. 
%
% -> f is defined as @ function, so that given a matrix M of points (each column being a point), f(M) gives a vector of values
%
% res = QUADRATURE_ON_SPARSE_GRID(f,S,map)
%
% -> map modifies points in (-1,1)^N to (a,b)^N and it's defined as @ function , that transforms a matrix (each column
%       is a point) in a matrix
%
% res = QUADRATURE_ON_SPARSE_GRID(f,S,map,weights_fact)
%
% -> weights_fact is a factor to modify the weights of the grid. If no map is needed, set map=[]


if exist('map','var') && ~isempty(map) % this works, && is shortcircuited
    knots=map([S.knots]);
else
    knots=[S.knots];
end

if exist('weights_fact','var') && ~isempty(weights_fact)
    weights=weights_fact*[S.weights]';
else
    weights=[S.weights]';
end

res=f(knots)*weights;