function [iseq,whatfield] = isequal_tensor_grids(T1,T2,Tol)


% ISEQUAL_TENSOR_GRID compares two tensor grids field by field, checking that the two grids are identical
%       and accounting for numerical tolerance when comparing knots and weights.
% 
%
% [ISEQ,WHATFIELD] = ISEQUAL_TENSOR_GRIDS(T1,T2) checks if sparse grids S1 and S2 are equal in the 'grid-by-grid' sense
%       The default tolerance 1e-14 is used when comparing points and weights.
%       The outputs are:
%       ISEQ = true / false
%       WHATFIELD = the name of the first field of the two tensor grids that are found different, or 'length' if the
%                   two sparse grids have a different number of tensor grids                    
%
%
% [ISEQ,WHATFIELD] = ISEQUAL_TENSOR_GRIDS(T1,T2,TOL) specifies the tolerance to compare knots and weights


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if nargin == 2
    Tol = 1e-14;
end


if is_tensor_grid(T1)==0 || is_tensor_grid(T2)==0 % (1 is for tensor, -1 is for tensor as parts of sparse grids)
    error('SparseGKit:WrongInput','T1 or T2 not tensor grids in isequal_tensor_grids')
end

if ~isequal(T1.size,T2.size) 
    iseq = false;
    whatfield = 'size';
    return
end


if ~isequal(T1.m, T2.m)
    iseq = false;
    whatfield = 'm';
    return
end


if max(abs((T1.weights - T2.weights))) > Tol
    iseq = false;
    whatfield = 'weights';
    return    
end    

N = length(T1.knots_per_dim);
N2 = length(T2.knots_per_dim);

if N~=N2
    iseq = false;
    whatfield = 'knots_per_dim';
    return
end

for n = 1:N
    if max(abs(T1.knots_per_dim{n} - T2.knots_per_dim{n})) > Tol       
        iseq = false;
        whatfield = 'knots_per_dim';
        return
    end
end
 
if max(max(abs((T1.knots - T2.knots)))) > Tol
    iseq = false;
    whatfield = 'knots';
    return    
end    


iseq = true;
whatfield = '';