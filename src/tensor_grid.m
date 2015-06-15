function S = tensor_grid(N,m,knots)

%TENSOR_GRID generates a tensor grid and computes the corresponding weights
%
%S = TENSOR_GRID(N,M,KNOTS) creates a tensor grid in N dimensions with M=[m1,m2,...,m_N] points 
%       in each direction. KNOTS is either a cell array containing the functions to be used 
%       to generate the knots in each direction, i.e.         
%
%            KNOTS={@knots_function1, @knots_function2, ... }
%
%       or a single function, to be used to generate the 1D knots in every direction, i.e.
%
%            KNOTS=@knots_function1
%
%       In both cases, the header of knots_function is [x,w]=knots_function(m)
%
%       The output S is a structure containing the information on the tensor grid:
%           S.knots: vector containing the tensor grid knots
%           S.weights: vector containing the corresponding weights
%           S.size: size of the tensor grid = prod(m)



%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% if knots is a simple function, we replicate it in a cell
if isa(knots,'function_handle')
    fknots=knots;
    knots=cell(1,N);
    for i=1:N
        knots{i}=fknots;
    end
end

sz=prod(m);

S.knots=zeros(N,sz);
S.weights=ones(1,sz);

% generate the pattern that will be used for knots and weights matrices, e.g.
% 
% pattern = [1 1 1 1 2 2 2 2;
%            1 1 2 2 1 1 2 2;
%            1 2 1 2 1 2 1 2]
%
% meaning "first node d-dim uses node 1 in direction 1, 2 and 3, second d-dim  node uses node 1 in 
% direction 1 and 2 and node 2 in direction 3 ...

pattern=generate_pattern(m);

for n=1:N
    [xx,ww] = knots{n}(m(n));
    S.knots(n,:) = xx(pattern(n,:));
    S.weights = S.weights.*ww(pattern(n,:));
end
S.size = sz;