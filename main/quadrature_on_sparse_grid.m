function [res,evals] = quadrature_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral,tol)

% QUADRATURE_ON_SPARSE_GRID uses a sparse grid to compute the integral of a function.
% It provides evaluation recycling and parallel toolbox support. This function behaves as
% EVALUATE_ON_SPARSE_GRID, except that it return the value of the approximated integral 
% of the function. See EVALUATE_ON_SPARSE_GRID for more information on inputs. Possible calls:
%
% res = QUADRATURE_ON_SPARSE_GRID(F,SR) where F is a function handle
%
% res = QUADRATURE_ON_SPARSE_GRID(F_VALS,SR) where F_VALS is a vector containing the evaluations of F over SR 
%       F can also be a matrix containing in each row the evaluation of a different function
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD) 
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,[],[],[]) 
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,[],SR_OLD) 
%
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD,PARAL) 
%
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD,PARAL,TOL)
%
%
% [res,evals] = QUADRATURE_ON_SPARSE_GRID(...) returns the evaluations of the function F 
%       over the points of the sparse grid S


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end


switch nargin
    
    case 1
        error('SparseGKit:WrongInput','not enough input arguments')
        
    case 2
        % res = QUADRATURE_ON_SPARSE_GRID(f,S), S being a reduced sparse grid. 
        if isa(f,'function_handle') && isreduced(S) 
            evals = evaluate_on_sparse_grid(f,S);
            res = evals*S.weights';
            
        elseif isnumeric(f) && isreduced(S) 
            
            res = f*S.weights';
            
        else
            error('SparseGKit:WrongInput','when quadrature_on_sparse_grid is called with two inputs, the second one must be a reduced sparse grid')
        end        
        return
        
    case {3,4,5}
        errmsg = ['QUADRATURE_ON_SPARSE_GRID does not accept ',num2str(nargin),' inputs.'];
        error('SparseGKit:WrongInput',errmsg)
                
    case 6
        evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old);
        
    case 7
        evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral);
        
    case 8
        evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral,tol);
end

res = evals*Sr.weights';