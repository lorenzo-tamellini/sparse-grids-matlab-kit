function [f_eval,new_points,tocomp_list,discard_points,discard_list] = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral,tol)

%EVALUATE_ON_SPARSE_GRID evaluates a function on a sparse grid, possibly recycling previous calls.
%           Several input combinations are possible, every input from the 3rd on can be set to []
% 
% F_EVAL = EVALUATE_ON_SPARSE_GRID(F,SR) evaluates a function F on a sparse grid SR, without recycling. 
%           SR can be a reduced sparse grid or a matrix with points stored as columns (i.e., like points
%           are stored in Sr.knots, Sr being a reduced sparse grid).
%           F is a function that takes as input a column vector point and returns either a scalar or a column vector. 
%           F will be evaluated one point at a time so there's no need for F to accept as input matrices as well.
%
%           F_EVAL is a matrix storing the evaluations of F on SR, each evaluation being stored as a column 
%           vector 
%
%
% F_EVAL = EVALUATE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD) recycles available evaluations of F on a different
%           sparse grid, stored respectively in EVALS_OLD and SR_OLD. EVALS_OLD is a matrix storing the 
%           evaluations of F on SR_OLD, each evaluation being stored as a column vector 
%           (i.e. EVALS_OLD will be typically a row vector or a matrix with nb.columns = nb. sparse grid
%           points). 
%
%           Two kinds of inputs are possible:
%
%           -> S and S_OLD are sparse grids, and SR, SR_OLD their reduced version, i.e.
%               SR=REDUCE_SPARSE_GRID(S) and analogous for SR_OLD
% 
%           -> SR and SR_OLD are just "list of points", i.e. matrices with points stored as columns,
%               S=[] and S_OLD=[]. This kind of call will be *much* slower than the other one though, especially for large N
%
%
% [F_EVAL,NEW_POINTS,IDX_NEW] = EVALUATE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD) also returns NEW_POINTS, the list of points 
%           where f has been evaluated (i.e. the new points w.r.t. the previous grid), and IDX_NEW, that contains
%           the position of NEW_POINTS in SR.KNOTS, i.e. SR.KNOTS(:,IDX_NEW) = NEW_POINTS
%
% [F_EVAL,NEW_POINTS,IDX_NEW,DISCARD_POINTS,IDX_OLD] = EVALUATE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD) also returns
%           DISCARD_POINTS, the list of points of SR_OLD that have been discarded (may happen with non-nested
%           grids) and IDX_OLD, index vector s.t. SR_OLD.KNOTS(:,IDX_OLD) = DISCARD_POINTS
%
%
% F_EVAL = EVALUATE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD,PARAL) uses the matlab parallel toolbox to speed 
%           up the computation:
%               PARAL = NaN means no parallel (default)
%               PARAL = some number X means "use parallel if more than X evals are needed". 
%                   This is useful for fast evals, in which case parallel may actually be inefficient 
%                    (due to communication time)
%           Important:  only function evaluations are processed in parallel; the preliminary
%           analysis (i.e. looking for recyclable evaluations) is performed in serial
%           Important: EVALUATE_ON_SPARSE_GRID does not switch on a matlabpool session. However, an error
%           is thrown if no matlabpool session is detected
%
%
% F_EVAL = EVALUATE_ON_SPARSE_GRID(F,S,SR,[],[],[],PARAL) uses the matlab parallel toolbox without recycling  
%
%
% F_EVAL = EVALUATE_ON_SPARSE_GRID(F,S,SR,EVALS_OLD,S_OLD,SR_OLD,PARAL,TOL) specifies the tolerance to be used when 
%            testing whether two points are equal or not (default 1e-14)
%
%
% F_EVAL = EVALUATE_ON_SPARSE_GRID(F,S,SR,[],[],[],PARAL,TOL) uses the matlab parallel toolbox without recycling  


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% declare a global variable controlling verbosity
global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end


% safety check and input handling
% ------------------------------------


% switch on nargin to decide whether to use simple_evaluate or go with the recycling. Cases of the switch
% range from 2 to 7 (8 is full input). Note that case 8 necessarely means that we are going for the recyle (if not, why
% specifying tol?)

switch nargin

    case 1
        error('not enough input arguments')

    case 2
        % evaluate_on_sparse_grid(f,Sr). 
        % The second input is S though, no we need to change its name
        if ~isreduced(S)
            Sr=as_reduced(S);
        else
            Sr=S;
        end
        f_eval = simple_evaluate(f,Sr);
        new_points = Sr.knots;
        tocomp_list = 1:length(Sr.weights);
        discard_points=[];
        discard_list=[];
        return

    case {3,4,5}
        errmsg = ['EVALUATE_ON_SPARSE_GRID does not accept ',num2str(nargin),' inputs. ' ... 
            'Observe that EVALUATE_ON_SPARSE_GRID has been changed after release 15.8 and ' ...
            'now takes as input also the non-reduced versions of the sparse grids on which ' ...
            'one wants to evaluate the function. This allows for significant savings in computational '...
            'time, especially for N large. See EVALUATE_ON_SPARSE_GRID for more information. '...
            'As a quick fix, you can use EVALUATE_ON_SPARSE_GRID_LEGACY '...
            'which is the old version of EVALUATE_ON_SPARSE_GRID, see help EVALUATE_ON_SPARSE_GRID_LEGACY. '...
            'This function is however deprecated and will disappear from future relesases of the Sparse Grid Matlab Kit.'];
        error(errmsg)

    case 6
        % evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old)
        % or
        % evaluate_on_sparse_grid(f,S,Sr,[],[],[])
        
        % in the previous version, the 6-input cases were
        %
        % evaluate_on_sparse_grid(f,Sr,[],[],paral,tol)
        % or
        % evaluate_on_sparse_grid(f,Sr,evals_old,Sr_old,paral,tol)
        % 
        % raise an error in these cases, which can be reckognized because the second input is 
        % not a sparse grid but a reduced sparse grid already
        if ~issmolyak(S)
            errmsg = ['The second input of EVALUATE_ON_SPARSE_GRID must be a non-reduced sparse grid if the function is called with 6 inputs. ' ...
                'Observe that EVALUATE_ON_SPARSE_GRID has been changed after release 15.8 and ' ...
                'now takes as input also the non-reduced versions of the sparse grids on which ' ...
                'one wants to evaluate the function. This allows for significant savings in computational '...
                'time, especially for N large. See EVALUATE_ON_SPARSE_GRID for more information. '...
                'As a quick fix, you can use EVALUATE_ON_SPARSE_GRID_LEGACY '...
                'which is the old version of EVALUATE_ON_SPARSE_GRID, see help EVALUATE_ON_SPARSE_GRID_LEGACY. '...
                'This function is however deprecated and will disappear from future relesases of the Sparse Grid Matlab Kit.'];
            error(errmsg)
        end
        
        if ~isreduced(Sr)
            Sr=as_reduced(Sr);
        end
        if isempty(evals_old) && isempty(Sr_old) && isempty(S_old)
            f_eval = simple_evaluate(f,Sr);
            new_points = Sr.knots;
            tocomp_list = 1:length(Sr.weights);
            discard_points=[];
            discard_list=[];
            return
        end
        if ~isreduced(Sr_old)
            error('SR_OLD must be a reduced sparse grid')
        end
        paral=NaN;
        tol=1e-14;

    case 7
        % evaluate_on_sparse_grid(f,S,Sr,[],[],[],paral)
        % or
        % evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old,paral)
        if ~isreduced(Sr)
            Sr=as_reduced(Sr);
        end
        if isempty(evals_old) && isempty(Sr_old)
            f_eval = simple_evaluate(f,Sr,paral);
            new_points = Sr.knots;
            tocomp_list = 1:length(Sr.weights);
            discard_points=[];
            discard_list=[];
            return
        end
        if ~isreduced(Sr_old)
            as_reduced(Sr_old);
        end
        tol=1e-14;

    case 8
        % evaluate_on_sparse_grid(f,S,Sr,[],[],[],paral,tol)
        % or
        % evaluate_on_sparse_grid(f,Sr,evals_old,Sr_old,paral,tol)
        if ~isreduced(Sr)
            Sr=as_reduced(Sr);
        end
        if isempty(evals_old) && isempty(Sr_old)
            f_eval = simple_evaluate(f,Sr,paral);
            new_points = Sr.knots;
            tocomp_list = 1:length(Sr.weights);            
            discard_points=[];
            discard_list=[];
            return
        end
        if ~isreduced(Sr_old)
            as_reduced(Sr_old);
        end        
end


% look for points that need to be evaluated
% ------------------------------------


pts_list = Sr.knots';
pts_list_old = Sr_old.knots';

% we store the needed information in three vectors:
%
% -> tocomp_list contains the indices of the points where f needs to be evaluated
% -> recycle_list contains the indices of the points that have been already evaluated
% -> recycle_list_old contains their position in the old sparse grid

if isempty(S) && isempty(S_old) % we go for the slow code
    [tocomp_list,recycle_list,recycle_list_old,discard_list] = lookup_merge_and_diff(pts_list,pts_list_old,tol);
elseif ~isempty(S) && ~isempty(S_old)
    [tocomp_list,recycle_list,recycle_list_old,discard_list] = compare_sparse_grids(S,Sr,S_old,Sr_old,tol);
else
    error('unreckongnized set of inputs')
end
    
new_points=pts_list(tocomp_list,:)';
discard_points = pts_list_old(discard_list,:)';

% do the evaluation 
% ------------------------------------


N=size(pts_list,1);
s = size(evals_old,1);

f_eval=zeros(s,N);

n = size(tocomp_list,1);
evals_new=zeros(s,n);

if ~isempty(tocomp_list)
    if n>paral % if no parallel this one becomes n>NaN, which is false for any n
        if MATLAB_SPARSE_KIT_VERBOSE,   
            disp('using parallel')
        end
        if ~check_if_parallel_on()
            error('no open matlabpool session detected')
        end
        parfor i=1:n
            % suppress the "variable is indexed but not sliced" warning, which cannot be circumvented in this case
            evals_new(:,i)=f(pts_list(tocomp_list(i),:)'); %#ok<PFBNS>
        end
    else
        % either no parallel available or no parallel wanted, so we just go with a for
        if MATLAB_SPARSE_KIT_VERBOSE
            disp('using serial')
        end
        for i=1:n
            evals_new(:,i)=f(pts_list(tocomp_list(i),:)');
        end
    end
end

% the two parts together 
% ------------------------------------

f_eval(:,tocomp_list)= evals_new;
f_eval(:,recycle_list)= evals_old(:,recycle_list_old);


end





%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------






function output = simple_evaluate(f,Sr,paral)


% function output = simple_evaluate(f,Sr,paral)
%
% does not exploit the previous sparse grids evaluations

% declare a global variable controlling verbosity

global MATLAB_SPARSE_KIT_VERBOSE
if isempty(MATLAB_SPARSE_KIT_VERBOSE)
    MATLAB_SPARSE_KIT_VERBOSE=1;
end


if nargin==2,
    paral=NaN;
end
    
n=size(Sr.knots,2);

% probe f to see the size of its output and 
probe=f(Sr.knots(:,1));
output=zeros(length(probe),n);
output(:,1)=probe;

if n>paral % if no parallel this one becomes n>NaN, which is false for any n
    if MATLAB_SPARSE_KIT_VERBOSE
        disp('using parallel')
    end
    if ~ check_if_parallel_on()
        error('no open matlabpool session detected')
    end
    parfor i=2:n
        % if ~mod(i,100), disp(i), end        
        output(:,i)=f(Sr.knots(:,i)); %#ok<PFBNS> % suppress the "variable is indexed but not sliced" warning, which cannot be circumvented in this case 
    end    
else
    % either no parallel available or no parallel wanted, so we just go with a for
    if MATLAB_SPARSE_KIT_VERBOSE
        disp('using serial')
    end
    for i=2:n
        % if ~mod(i,100), disp(i), end
        output(:,i)=f(Sr.knots(:,i)); 
    end    
end
    
    
end

