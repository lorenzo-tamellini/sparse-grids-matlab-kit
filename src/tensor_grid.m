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




% if knots is a simple function, we replicate it in a cell
if isa(knots,'function_handle')
    fknots=knots;
    knots=cell(1,N);
    for i=1:N
        knots{i}=fknots;
    end
end

[S.knots,S.weights] = knots{1}(m(1));   
for i=2:N
 [Xi,Wi] = knots{i}(m(i));
 S.knots = combvec(S.knots,Xi);         
 S.weights = combvec(S.weights,Wi);
 S.weights = prod(S.weights);
end

S.size = prod(m);
