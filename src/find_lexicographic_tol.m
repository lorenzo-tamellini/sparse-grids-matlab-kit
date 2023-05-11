function [found,pos,iter] = find_lexicographic_tol(lookfor,I,TOL)

% FIND_LEXICOGRAPHIC finds specific rows of a matrix that is sorted lexicographically up to TOL. 
%   A tolerance TOL to test equality between LOOKFOR and rows of I, and also to test whether
%   LOOKFOR and rows of LOOKFOR  are lexico-sorted. 
%
% [FOUND,POS] = FIND_LEXICOGRAPHIC_TOL(LOOKFOR,I) looks for the row vector LOOKFOR among
%               the rows of I. If found, FOUND=TRUE and POS is the rownumber of LOOKFOR,
%               i.e. LOOKFOR == I(POS,:) up to tolerance 1e-14. If not found, FOUND=FALSE and POS=[]. 
%
%               The function **does not** performs a preliminar check whether I is actually lexicographically sorted.
%   
% [FOUND,POS] = FIND_LEXICOGRAPHIC_TOL(LOOKFOR,I,TOL) sets internal tolerance to TOL
%
%
% Es:
% I = multiidx_box_set([4 2],1); %I contains [4 2]
% noise = 1e-13*randn(size(I));
% Inoisy = I + noise;
% 
% find_lexicographic_tol([4 2],Inoisy,1e-12) % -> true
% find_lexicographic_tol([4 2],Inoisy,1e-14) % -> false




%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% fixing optional inputs
if nargin == 2
    TOL = 1e-14;
end


nb_idx = size(I,1);

% we exploit the fact that I is sorted lexicographically and we proceed by binary search
% which guarantees cost ~ log(nb_idx)
%
% Basically, we start from the middle row, compare it with the index to be found, and
% if our index is larger we make a step in the increasing direction (i.e. we look in the upper half
% of the sorting), and viceversa. Of course, the step halves at each iteration: 
% therefore we necessarily terminate in ceil(log2(nb_idx)) steps at most

% the position to compare against -- if found, this is the position to be returned
idx = ceil(nb_idx/2);
% the associated vector 
jj = I(idx,:); 

found = isequal_tol(jj,lookfor,TOL);

iter=1;
itermax = ceil(log2(nb_idx));

while ~found && iter <= itermax
    if islexico_tol(jj,lookfor,TOL) % look in the upper half, unless we are already at the end
        if idx < nb_idx
            idx = min(idx+ceil(nb_idx/2^iter),nb_idx);
            jj = I(idx,:); % row entry of C2
            found = isequal_tol(lookfor,jj,TOL);
            iter=iter+1;
        else
            break
        end
    else % look in the bottom half, unless we are already at the beginning
        if idx > 1
            idx = max(idx-ceil(nb_idx/2^iter),1);
            jj = I(idx,:); % row entry of C2
            found = isequal_tol(lookfor,jj,TOL);
            iter=iter+1;
        else
            break
        end
    end
end

pos = [];
if found 
    pos = idx;
end

end




