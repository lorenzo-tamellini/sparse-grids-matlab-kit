function [iseq,whatfield] = isequal_sparse_grids(varargin)

% ISEQUAL_SPARSE_GRID compares two sparse grids, checking that the two grids are identical tensor grid by tensor grid.
%      In practice, a version of the standard ISEQUAL function for vectors of structs, where however
%      we account for numerical tolerance when comparing knots and weights.
% 
%
% [ISEQ,WHATFIELD] = ISEQUAL_SPARSE_GRIDS(S1,S2) checks if sparse grids S1 and S2 are equal in the 'grid-by-grid' sense
%       The default tolerance 1e-14 is used when comparing points and weights.
%       The outputs are:
%       ISEQ = true / false
%       WHATFIELD = the name of the first field of the two tensor grids that are found different, or 'length' if the
%                   two sparse grids have a different number of tensor grids                    
%
%
% [ISEQ,WHATFIELD] = ISEQUAL_SPARSE_GRIDS(S1,S2,TOL) specifies the tolerance to compare knots and weights


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



switch nargin

    case 2 % iseq = isequal_sparse_grids(S1,S2)

        S1 = varargin{1};
        S2 = varargin{2};
        Tol = 1e-14;

        if ~is_sparse_grid(S1) || ~is_sparse_grid(S2)
            error('SparseGKit:WrongInput','S1 or S2 not sparse grids in isequal_sparse_grids')
        end
        
        
    case 3 % iseq = isequal_sparse_grids(S1,S2,Tol)
        
        S1 = varargin{1};
        S2 = varargin{2};
        Tol = varargin{3};

        if ~is_sparse_grid(S1) || ~is_sparse_grid(S2)
            error('SparseGKit:WrongInput','S1 or S2 not sparse grids in isequal_sparse_grids')
        end        
        
    otherwise
        error('SparseGKit:WrongInput','unknown inputs for isequal_sparse_grids')
end


    
L1 = length(S1);
L2 = length(S2);

if L1 ~= L2
    iseq = false;
    whatfield = 'length';
    return
end

for l =  1:L1
    
    [iseq,whatfield] = isequal_tensor_grids(S1(l),S2(l),Tol);
    
    if ~iseq
        return
    end
    
    if ~isequal(S1(l).coeff,S2(l).coeff)
        iseq = false;
        whatfield = 'coeff';
        return
    end
    
    if ~isequal(S1(l).idx,S2(l).idx)
        iseq = false;
        whatfield = 'idx';
        return
    end
    
end


