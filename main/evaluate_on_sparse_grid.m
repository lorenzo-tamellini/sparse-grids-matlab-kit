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







%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------






function [tocomp_list,recycle_list,recycle_list_old] = lookup_merge_and_diff(pts_list,pts_list_old,Tol)

% [tocomp_list,recycle_list,recycle_list_old] = lookup_merge_and_diff(pts_list,pts_list_old,tol)
%
% looks for points of pts_list in pts_list_old using the same algorithm as reduce_sparse_grid. tol is
% the tolerance for two points to be considered equal

if nargin==2
    Tol=1e-14;
end

N=size(pts_list,1);
N_old=size(pts_list_old,1);

% the list of indices of points to be evaluated. Init to its max length, will be cut after the search is over
tocomp_list=zeros(N,1);
% if grid are not nested, not all of the old grid will be recycled. Thus,
% we also need 2 list of indices of points to be recycled, storing the
% positions  in the old grid and in the new one
recycle_list_old=zeros(N_old,1);
recycle_list=zeros(N,1);



% first, I merge the two lists 

Merged = [pts_list_old; pts_list];

% and order the rows of Merged in lexicographic order, obtaining Sorted. If I use mysortrows then two rows like
%
% [a b c d]
% [a-t b c+t d]
% 
% are considered equal (if t < Tol ) and therefore placed one after the other
%
% sorter is an index vector that maps Merged into Sorted, i.e. Merged(sorter,:)==Sorted

[Sorted,sorter] = mysortrows(Merged,Tol/size(Merged,2));

% I also need to remember which points come from pts_list and which from
% pts_list_old, and from which original position. 
% Thus I create a flag vector [-1 -2 ... -N_old 1 2 .. N], i.e. positive
% flags identify the new grid and negative ones the old grid, then I sort
% it according to sorter

flags=[-1:-1:-N_old, 1:N];
flags_sorted=flags(sorter);


% next I take the difference of two consecutive rows. if the difference is small, then the rows are the same, i.e. the knot is the same
dSorted=diff(Sorted,1,1);

% I measure the difference with infty norm instead of L2 norm:
% i take  the maximum component of each row (2 means "operate on columns"):
%       max(abs(dSorted),[],2)    
% then I want to see which ones have this max bigger than Tol
%       diff_eq=max(abs(dSorted),[],2)>Tol
% this command returns a vector of True and false
%       diff_eq=[1 1 0 1 1]
% this means that the 2nd point is different from the 1st, the 3rd from the 2nd, but
% the 4th is equal to the 3rd ( diff(3)=v(4)-v(3) ), hence in common between the
% grids.

diff_eq=max(abs(dSorted),[],2)>Tol;

% now I scroll diff_eq and sort out everything according to these rules:
%
% --> if diff_eq(k)==0, 
%
% in this case either Sorted(k+1) is in the new grid and Sorted(k)
% is in the old grid and both are equal or viceversa. 
%
% Therefore, the point in the new grid goes into recycle_list and the
% old one in recycle_list_old. 
%
% Then I can skip the following (since it's equal and I have sorted it
% already)
%
% --> else, diff_eq(k)==1. 
%
% In this case, 4 cases are possible but actually only two matters
%
% -----> both Sorted(k) and Sorted(k+1) comes from the old_grid. 
%   then, Sorted(k) is to discard
%
% -----> both Sorted(k) and Sorted(k+1) comes from the new_grid. 
%   then, Sorted(k) is to compute
%
% -----> Sorted(k) is new, Sorted(k+1) is old. 
%   then, Sorted(k) is to compute
%
% -----> Sorted(k) is old, Sorted(k+1) is new. 
%   then, Sorted(k) is to discard



i=0; % scrolls recycle lists
j=0; % scroll compute_list
k=1; % scrolls diff_eq
L=length(diff_eq);
discard=0;

while k<=L
    
    if diff_eq(k) % short-hand for diff_eq(k)==1, i.e. two consecutive entries are different
        
       % compute or discard        
       if flags_sorted(k)>0
           j=j+1; 
           tocomp_list(j)=flags_sorted(k);
       else
           discard=discard+1;
       end
       
       % then move to the following
       k=k+1;
       
    else % in this case diff_eq(k)==0, i.e. two consecutive entries are equal
         
        % recycling case
        i=i+1;
        
        if flags_sorted(k)>0    
            recycle_list(i)=flags_sorted(k);
            recycle_list_old(i)=-flags_sorted(k+1);
        else
            recycle_list_old(i)=-flags_sorted(k);
            recycle_list(i)=flags_sorted(k+1);
        end
        
        % then I can skip the k+1 because I have already sorted it
        k=k+2;
    end
    
end

% need to handle the case k==L, since diff is 1 element shorter than sorted. Note
% that 
% --> the node in Sorted(L,:) has been already taken care of inside the while loop
% --> we only need to do something if the diff_eq(L)==1. Indeed, if 
% diff_eq(L)==0, then Sorted(L+1,:) is already taken care of (and in this
% case the final value of k is L+2)

if diff_eq(L) % short-hand for diff_eq(L)==1,
    if flags_sorted(L+1)>0 % then it's a new point and has to be computed
        j=j+1;
        tocomp_list(j)=flags_sorted(L+1);
    else % then it's an old point and has to be discarded
        discard=discard+1;
    end
end

% show some statistics
display(strcat('new evaluation needed:',num2str(j),' recycled evaluations:',num2str(i),' discarded evaluations:',num2str(discard)))



% remove the extra entries of tocomp_list and recycle_lists
if tocomp_list(j+1)~=0,error('tocomp_list(j+1)~=0'),end
tocomp_list(j+1:end)=[];
if recycle_list(i+1)~=0,error('recycle_list(i+1)~=0'),end
recycle_list(i+1:end)=[];
if i~=length(recycle_list_old) % in this case we haven't exhausted the old_grid (non-nested case)
    if recycle_list_old(i+1)~=0,error('recycle_list_old(j+1)~=0'),end
    recycle_list_old(i+1:end)=[];
end

% safety checks
if length(recycle_list)~=length(recycle_list_old),
    error('length(recycle_list)~=length(recycle_list_old)'),
end
if ~isempty( setxor( [tocomp_list; recycle_list], 1:N ) ),
    error('~isempty(setxor([tocomp_list recycle_list],1:N))'),
end


end