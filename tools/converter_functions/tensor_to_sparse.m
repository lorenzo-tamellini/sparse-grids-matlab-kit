function S = tensor_to_sparse(T,idx)

% TENSOR_TO_SPARSE converts a tensor grid into a sparse grid structure of the same type of CREATE_SPARSE_GRID, by adding the missing fields
% 
% S = TENSOR_TO_SPARSE(T) copies the fields of T in S and adds the following
%           S.coeff = 1
%           S.idx   = vector with idx(i) = length(T.knots_per_dim{i})
%
% S = TENSOR_TO_SPARSE(T,idx) copies the fields of T in S and adds the following
%           S.coeff = 1
%           S.idx   = idx

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------


S = T;  
S.coeff = 1; 
if nargin == 2
    S.idx = idx;  
else
    N = length(S.knots_per_dim);
    idx = zeros(1,N);
    for n = 1:N
        idx(n) = length(S.knots_per_dim{n});
    end
    S.idx = idx;
end
