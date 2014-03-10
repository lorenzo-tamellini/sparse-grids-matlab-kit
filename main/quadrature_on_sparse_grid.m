function res = quadrature_on_sparse_grid(f,S,evals_old,S_old,paral,tol)

% QUADRATURE_ON_SPARSE_GRID uses a sparse grid to compute the integral of a function.
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S) computes the integral of F using the sparse grid S.
%       S is a reduced sparse grid; F is a function that takes as input a column vector 
%       (i.e. a sparse grid point) and returns either a scalar or a
%       column vector. F will be evaluated one point at a time so there's no need for F to accept as input 
%       matrices as well.
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,EVALS_OLD,S_OLD) recycles available evaluations 
%       of F on a different sparse grid (SR_OLD). EVALS_OLD is a matrix storing the evaluations 
%       of F on SR_OLD, where each evaluation stored as a column vector; SR_OLD is a reduced sparse_grid 
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,EVALS_OLD,SR_OLD,PARAL) specifies the parameter
%       controlling the use of the parallel computing toolbox, see EVALUATE_ON_SPARSE_GRID. 
%
% res = QUADRATURE_ON_SPARSE_GRID(F,S,EVALS_OLD,S_OLD,PARAL,TOL), specifies the tolerance
%       used when checking equality between points of S and S_OLD. 
       

switch nargin
    case 1
        error('not enough input arguments')
    case 2
        % res = QUADRATURE_ON_SPARSE_GRID(f,S)
        evals = evaluate_on_sparse_grid(f,S);
    case 3
        % this was once this call, no longer valid
        % res = QUADRATURE_ON_SPARSE_GRID(f,S,map)
        errmsg=['the latest version of QUADRATURE_ON_SPARSE_GRID does not accept 3 input arguments. '...
                'In particular, the call QUADRATURE_ON_SPARSE_GRID(f,S,map) is no longer valid. '...
                'To fix this, either regenerate your sparse grid using MAP as input (help SMOLYAK_GRID) '...
                'or modify the field KNOTS of your reduced sparse grid applying MAP '...
                '(e.g.  Sr.knots = map(Sr.knots) ). '...
                'This message will disappear in future releases of SPARSE_GRID_MATLAB_KIT'];
        error(errmsg);
    case 4
        % this is either
        % % res = QUADRATURE_ON_SPARSE_GRID(f,S,map,weights_fact)
        % which is no longer valid, or
        % res = QUADRATURE_ON_SPARSE_GRID(f,S,evals_old,S_old)
        % or
        % res = QUADRATURE_ON_SPARSE_GRID(f,S,[],[])
        if isa(evals_old,'function_handle') || isscalar(S_old)
            errmsg=['The call QUADRATURE_ON_SPARSE_GRID(f,S,map,weights_fact) is no longer valid. '...
                'To fix this, either regenerate your sparse grid using MAP  and WEIGHTS_FACT as input (help SMOLYAK_GRID) '...
                'or modify the fields KNOTS and WEIGHTS of your reduced sparse grid applying MAP and WEIGHTS_FACT'...
                '(e.g.  Sr.knots = map(Sr.knots),  Sr.weights=weights_fact*Sr.weights). '...
                'This message will disappear in future releases of SPARSE_GRID_MATLAB_KIT'];
            error(errmsg);
        end
        
        evals = evaluate_on_sparse_grid(f,S,evals_old,S_old);
        
    case 5
        evals = evaluate_on_sparse_grid(f,S,evals_old,S_old,paral);
    case 6
        evals = evaluate_on_sparse_grid(f,S,evals_old,S_old,paral,tol);

end

% o is the number of outputs
s=size(evals,1);

res=zeros(s,1);
for d=1:s
    res(s)=sum(dot(evals(s,:),S.weights));
end
