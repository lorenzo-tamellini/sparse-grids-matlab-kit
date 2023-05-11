function iss = issmolyak(S)

% iss = issmolyak(S)
% 
% DEPRECATED!! This function has been replaced by IS_SPARSE_GRID in release 23.5 and will disappear in future releases
%
% ISSMOLYAK(S) returns 1 if S is a sparse grid. A sparse grid is a vector of structs with fields 
% 'knots','weights','size','knots_per_dim','m','coeff','idx'.


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

error('SparseGKit:RenamedFunction', ['issmolyak has been replaced by is_sparse_grid in release 23.5.'...
    ' This message will disappear in future releases of the sparse-grid-matlab-kit.'])


if isstruct(S)
    iss=isempty(setxor(fieldnames(S),{'knots','weights','size','knots_per_dim','m','coeff','idx'})) && length(S)>=1;
else
    iss=0;
end

