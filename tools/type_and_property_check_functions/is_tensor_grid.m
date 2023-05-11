function ist = is_tensor_grid(S)

% IS_TENSOR_GRID(S) returns 1 or 0 or -1 
% -->  1 if returned if S is a tensor grid. A tensor grid is a struct with fields:
%           'knots','weights','size','knots_per_dim','m'
% --> -1 if returned if S is a tensor grid in a sparse grid struct. A tensor grid in a sparse grid struct is
%        a struct with two more fields:
%           'knots','weights','size','knots_per_dim','m','coeff','idx'
% --> 0 is returned in any other case


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

if isstruct(S)
    if isempty(setxor(fieldnames(S),{'knots','weights','size','knots_per_dim','m'})) && length(S)==1
        ist = 1;
    elseif isempty(setxor(fieldnames(S),{'knots','weights','size','knots_per_dim','m','coeff','idx'})) && length(S)==1
        ist = -1;
    else
        ist = 0;
    end
else
    ist=0;
end
