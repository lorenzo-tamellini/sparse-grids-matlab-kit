function [S,I,coeff] = create_sparse_grid_add_multiidx(new_idx,S_in,I_in,coeff_in,knots,lev2knots)

% CREATE_SPARSE_GRID_ADD_MULTIIDX produces a grid obtained by adding a single multi-idx to a previously existing grid.
%
% [S,I,COEFF] = CREATE_SPARSE_GRID_ADD_MULTIIDX(NEW_IDX,S_IN,I_IN,COEFF_IN,KNOTS,LEV2KNOTS) takes as inputs:
%       --> an index NEW_IDX. It must be admissible wrt to the index set I_IN described below, and this condition **won't** be checked
%           by the function
%       --> a sparse grid S_IN to which NEW_IDX should be added
%       --> the index set I_IN that was used to create S_IN, either implicitly (by defining the rule in CREATE_SPARSE_GRID)
%           or explicitely (by using CREATE_SPARSE_GRID_MULTIIDX). Note that this cannot be produced as the union of the indices S_IN.idx,
%           because S_IN contains only tensors whose coefficient in the combination technique is non-zero, while here we need them all.
%       --> COEFF_IN is the vector of coefficients of the combination technique obtained on I_IN. Note that this
%           vector of coefficients *MUST* include also zeros if a row of I_IN has coeff zero in the combination technique.
%           Therefore, COEFF_IN cannot be taken as [S_in.coeff]. Instead, use COEFF_IN = COMBINATION_TECHNIQUE(I_IN)           
%       --> KNOTS, LEV2KNOTS are the usual arguments to specify knots and lev2knots (see e.g. CREATE_SPARSE_GRID)
%
%       The function outputs:
%
%       --> S, the new sparse grid
%       --> I, the new multiidx set, i.e. I = sortrows([I_IN; NEW_IDX])
%       --> COEFF, the updated vector of coefficients of the combination technique


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


[I,sorter] = sortrows([I_in; new_idx]);

% add 0 to the new coeff_G in its right place (+1 will come later on)
coeff = [coeff_in 0];
coeff = coeff(sorter);


% first a check on C_old being sorted. Observe that the function sortrows used is very efficient
% so the cost of this preliminary analysis is negligible (e.g. it takes 0.02 sec to verify
% that a TD set with w=5 and N=30, i.e. ~340600 indices is sorted,  and only 0.00027 sec that
% the matrix A=randi(20,300000,30) is unsorted.

if ~issorted(I_in,'rows')
    error('SparseGKit:SetNotSorted','the multiindex set C_old is not sorted')
end

N=size(I_in,2);

% if knots and  lev2knots are simple function, we replicate them in a cell
if isa(knots,'function_handle')
    fknots=knots;
    knots=cell(1,N);
    for i=1:N
        knots{i}=fknots;
    end
end
if isa(lev2knots,'function_handle')
    f_lev2knots=lev2knots;
    lev2knots=cell(1,N);
    for i=1:N
        lev2knots{i}=f_lev2knots;
    end
end


%-----------------------------------------------
% generate contributions of jj to combitec. They are already lexicosorted, jj is the last row
Tensors_jj = delta_to_combitec(new_idx);


% look them up in G and for each of them add +-1. Because both Tensor and GG are sorted,
% restrict the search at each iteration

nb_tens = size(Tensors_jj,1);

from_where = 1;
for t = 1:nb_tens
    
    tt = Tensors_jj(t,:);
    [found,rel_pos] = find_lexicographic(tt,I(from_where:end,:),'nocheck'); % we know that G is sorted, so I can use the flag 'nocheck'
    
    pos = from_where + rel_pos - 1;
    
    if ~found
        error('I was supposed to find this tensor!')
    end
    
    coeff_tt = (-1)^sum(new_idx - tt);
    coeff(pos) = coeff(pos) + coeff_tt;
    
    from_where = pos + 1;
end


% now we can store only those grids who survived, i.e. coeff~=0
%------------------------------------------------------

nb_grids=sum(coeff~=0);
empty_cells=cell(1,nb_grids);
S=struct('knots',empty_cells,'weights',empty_cells,'size',empty_cells,'knots_per_dim',empty_cells,'m',empty_cells);
fieldnms = fieldnames(S)'; % I'll later need to loop over these fields - note the transpose, it has to be a row cell
coeff_condensed=zeros(1,nb_grids);
ss=1;

% for each nonzero coeff, generate the tensor grid and store it. If possible, recycle from S_old. Note that C_old
% stores all idx of S_old, even those with 0 coeff. So I need to extract the information on the tensors that I have
% in S_old as follows

if MATLAB_SPARSE_KIT_VERBOSE
    disp('build sparse grid with tensor grid recycling')
end

nb_S_old_grids = length(S_in);
C_old_nnz = reshape([S_in.idx],N,nb_S_old_grids)'; % [S_old.idx] puts all indices on a row. I reshape them and each N elements I get one column.
                                                    % Then I need to traspose the matrix



for j=1:length(coeff)
    if coeff(j)~=0
        i = I(j,:);
        [found,pos] = find_lexicographic(i,C_old_nnz,'nocheck'); % this exploits that C2 is lexicographic, so it's efficient (cost ~ log(nb_rows_C2))
        if found
            %disp('found')
            % Note that at this point elements of S are tensor grids while S_old is a sparse grid therefore it has additional fields
            % (coeff, idx). We thus need to copy field by field otherwise we'll have "assignment between dissimilar
            % structures" error. We use dynamic filed names to this end
            for fn = fieldnms
                S(ss).(fn{1}) = S_in(pos).(fn{1}); % note that each fn is a 1-element cell, so to access its value I need the notation fn{1}
            end
            % however we need to fix the weights. Indeed, they are stored in S_old as weights*coeff, so we need to reverse
            % that multiplication
            S(ss).weights = S(ss).weights/S_in(pos).coeff;
        else
            m =apply_lev2knots(i,lev2knots,N); % n. of points in each direction
            S(ss) = tensor_grid(N,m,knots);
        end
        S(ss).weights=S(ss).weights*coeff(j);
        coeff_condensed(ss)=coeff(j);
        ss=ss+1;
    end
end

    
% now store the coeff value. It has to be stored after the first loop, because tensor_grid returns a grid
% WITHOUT coeff field, and Matlab would throw an error (Subscripted assignment between dissimilar structures)

for ss=1:nb_grids
    S(ss).coeff=coeff_condensed(ss);
end

% similarly for the multiidx generating each tensor grid
ss=1;
for j=1:length(coeff)
    if coeff(j)~=0
        i = I(j,:);
        S(ss).idx = i;
        ss=ss+1;
    end
end


end % this end closes the function






function m = apply_lev2knots(i,lev2knots,N)
    
% N could be deduced by N but it's better passed as an input, to speed computation
% init m to zero vector
m=0*i;

% next, iterate on each direction 1,...,N. 

for n=1:N
    m(n) = lev2knots{n}(i(n));
end

end
