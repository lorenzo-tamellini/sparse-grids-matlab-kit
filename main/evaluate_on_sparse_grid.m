function output = evaluate_on_sparse_grid(f,Sr,evals_old,Sr_old,paral,tol)

%EVALUATE_ON_SPARSE_GRID evaluates a function on a sparse grid, possibly recycling previous calls
% 
% OUTPUT = EVALUATE_ON_SPARSE_GRID(F,SR) evaluates a function F on a sparse grid SR, without recycling. 
%           SR must be a reduced sparse grid. F is a function that takes as input a column vector point 
%           and returns either a scalar or a column vector. F will be evaluated one point at a time 
%           so there's no need for F to accept as input matrices as well.
%           OUTPUT is a matrix storing the evaluations of F on SR, each evaluation being stored as a column 
%           vector 
%
% OUTPUT = EVALUATE_ON_SPARSE_GRID(F,SR,EVALS_OLD,SR_OLD) recycles available evaluations of F on a different
%           sparse grid, stored respectively in EVALS_OLD and SR_OLD. EVALS_OLD is a matrix storing the 
%           evaluations of F on SR_OLD, each evaluation being stored as a column vector 
%           (i.e. EVALS_OLD will be typically a row vector or a matrix with nb.columns = nb. sparse grid points)
%
% OUTPUT = EVALUATE_ON_SPARSE_GRID(F,SR,EVALS_OLD,SR_OLD,PARAL) uses the matlab parallel toolbox to speed 
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
% OUTPUT = EVALUATE_ON_SPARSE_GRID(F,SR,[],[],PARAL) uses the matlab parallel toolbox without recycling  
%
% OUTPUT = EVALUATE_ON_SPARSE_GRID(F,SR,EVALS_OLD,SR_OLD,PARAL,TOL) specifies the tolerance to be used when 
%            testing whether two points are equal or not (default 1e-14)



% safety check and input handling
% ------------------------------------


if ~isreduced(Sr)
    error('SR must be a reduced sparse grid')
end


% switch on nargin to decide whether to use simple_evaluate or go with the recycling. Cases of the switch
% range from 2 to 5. Note that case 6 necessarely means that we are going for the recyle (if not, why
% specifying tol?)

switch nargin

    case 2
        % evaluate_on_sparse_grid(f,Sr)
        output = simple_evaluate(f,Sr);
        return

    case 4
        % evaluate_on_sparse_grid(f,Sr,evals_old,Sr_old)
        % or
        % evaluate_on_sparse_grid(f,Sr,[],[])
        if isempty(evals_old) && isempty(Sr_old)
            output = simple_evaluate(f,Sr);
            return
        end
        if ~isreduced(Sr_old)
            error('SR_OLD must be a reduced sparse grid')
        end
        paral=NaN;
        tol=1e-14;

    case 5
        % evaluate_on_sparse_grid(f,Sr,[],[],paral)
        % or
        % evaluate_on_sparse_grid(f,Sr,evals_old,Sr_old,paral)
        if isempty(evals_old) && isempty(Sr_old)
            output = simple_evaluate(f,Sr,paral);
            return
        end
        if ~isreduced(Sr_old)
            error('SR_OLD must be a reduced sparse grid')
        end
        tol=1e-14;

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

[tocomp_list,recycle_list,recycle_list_old] = lookup_merge_and_diff(pts_list,pts_list_old,tol);


% do the evaluation 
% ------------------------------------


N=size(pts_list,1);
s = size(evals_old,1);

output=zeros(s,N);

n = size(tocomp_list,1);
evals_new=zeros(s,n);

if ~isempty(tocomp_list)
    if n>paral % if no parallel this one becomes n>NaN, which is false for any n
        disp('solve with parallel')
        if ~matlabpool('size')
            error('no open matlabpool session detected')
        end
        parfor i=1:n
            % suppress the "variable is indexed but not sliced" warning, which cannot be circumvented in this case
            evals_new(:,i)=f(pts_list(tocomp_list(i),:)'); %#ok<PFBNS>
        end
    else
        % either no parallel available or no parallel wanted, so we just go with a for
        disp('solve with serial')
        for i=1:n
            evals_new(:,i)=f(pts_list(tocomp_list(i),:)');
        end
    end
end

% the two parts together 
% ------------------------------------


output(:,tocomp_list)= evals_new;
output(:,recycle_list)= evals_old(:,recycle_list_old);


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


if nargin==2,
    paral=NaN;
end

n=size(Sr.knots,2);

% probe f to see the size of its output and 
probe=f(Sr.knots(:,1));
output=zeros(length(probe),n);
output(:,1)=probe;

if n>paral % if no parallel this one becomes n>NaN, which is false for any n
    disp('solve with parallel')
    if ~matlabpool('size')
        error('no open matlabpool session detected')
    end
    parfor i=2:n
        % if ~mod(i,100), disp(i), end        
        output(:,i)=f(Sr.knots(:,i)); %#ok<PFBNS> % suppress the "variable is indexed but not sliced" warning, which cannot be circumvented in this case 
    end    
else
    % either no parallel available or no parallel wanted, so we just go with a for
    disp('solve with serial')
    for i=2:n
        % if ~mod(i,100), disp(i), end
        output(:,i)=f(Sr.knots(:,i)); 
    end    
end
    
    
end

