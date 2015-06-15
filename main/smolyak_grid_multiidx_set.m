function [S,C] = smolyak_grid_multiidx_set(C,knots,lev2knots,map,weights_coeff)

% SMOLYAK_GRID_MULTIIDX_SET produces a sparse grid starting from a multiindex-set rather than
% from a rule IDXSET(I) <= W.
%
% [S,C] = SMOLYAK_GRID_MULTIIIDX_SET(C,KNOTS,LEV2KNOTS) uses the multiindex set C. C must be
%       in lexicographic order and admissible. 
%
% [S,C] = SMOLYAK_GRID_MULTIIDX_SET(C,KNOTS,LEV2KNOTS,MAP,WEIGHTS_COEFF) can be used as an alternative
%       to generate a sparse grid on a hyper-rectangle.
%
%
% See also CHECK_SET_ADMISSIBILITY for admissibility check, and SMOLYAK_GRID for further information on
% KNOTS, LEV2KNOTS, MAP, WEIGHTS_COEFF and on the sparse grid data structure S


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


% first a check on C being sorted. Observe that the function sortrows used is very efficient
% so the cost of this preliminary analysis is negligible (e.g. it takes 0.02 sec to verify
% that a TD set with w=5 and N=30, i.e. ~340600 indices is sorted,  and only 0.00027 sec that
% the matrix A=randi(20,300000,30) is unsorted.

if ~issorted(C,'rows'),
    error('the multiindex set C is not sorted')
end

N=size(C,2);

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




N=size(C,2);

%-----------------------------------------------

% each multiindex of C is a delta grid. Given say [i1 i2 i3] i will get
% (i1 - i1-1) x (i2 - i2-1) x (i3 - i3-1)
% that is the grids associated to
% [i1 i2 i3],[i1-1 i2 13],[i1 i2-1 i3],[i1-1 i2-1 i3] and so on
%
% note that all of these multiindeces are also included in the initial set C
% (if [i1 i2 i3] ratisfies the rule, all other will, because their indeces are equal or lower)
%
% So the set of all possible grids (non delta grids) is the same set C, but some of them will cancel out
%
% C is partially ordered (lexicographic): [x x x x i], are listed increasing with i,
% [x x x j i] are listed increasing first with j and then with i ...
% Now take a row of C, c. Because of the ordering, if you take c as a grid index
% the same grid can appear again only from delta grids coming from rows following.
%
% Now we scroll these rows following (say c2) and compute how many times c will be generated, and with
% what sign. It has to be that d=c2-c has only 0 or 1
%
% if this is the case, the sign of this new appearence of c could be both + or - .
% To determine the sign, start with + and switch it every time a 1 appears in d
%
% d=[0 1 0 0 1] => sign= +
%
% the formula is (-1)^sum(d) ( d with even appearences of 1 gives a + , odd gives a -)
%
% and just sum all the coefficients to see if c will survive or cancel out
%
% so the algorithm works like this:

nn=size(C,1);
coeff=ones(1,nn); % initialize coefficients to 1: all c survive



% I can at least restrict the search to multiindices whose first component is c(i) + 2, so I define
[~,bookmarks]=unique(C(:,1),'first');
bk = [bookmarks(3:end)'-1 nn nn];
% i.e. those who begin with 1 end at bookmark(3)-1, those who begin with 2-1 end at bookmark(4) and so on,
% until there's no multiindex with c(i)+2
%
% an old piece of code I still want to have written somewhere
% range = i+ find( C(i+1:end,1)==cc(1)+2, 1 ) - 1;


for i=1:nn % scroll c
    cc = C(i,:);
    % recover the range in which we have to look. Observe that the first column of C contains necessarily 1,2,3 ...
    % so we can use them to access bk
    range=bk(cc(1));
    for j=i+1:range
        % scroll c2, the following rows
        d=C(j,:)-cc;
        if max(d)<=1 && min(d)>=0  % much faster to check then if only 0 and 1 appears. Also, && is short-circuited,
            % so if max(d)<=1 is false the other condition is not even checked
            coeff(i)=coeff(i) + (-1)^sum(d);
        end
    end
end



% now we can store only those grids who survived, i.e. coeff~=0
%------------------------------------------------------

nb_grids=sum(coeff~=0);
empty_cells=cell(1,nb_grids);
S=struct('knots',empty_cells,'weights',empty_cells,'size',empty_cells);
coeff_condensed=zeros(1,nb_grids);
ss=1;


for j=1:nn
    if coeff(j)~=0
        i = C(j,:);
        m =apply_lev2knots(i,lev2knots,N); % n. of points in each direction
        S(ss) = tensor_grid(N,m,knots);
        S(ss).weights=S(ss).weights*coeff(j);
        coeff_condensed(ss)=coeff(j);
        ss=ss+1;
    end
end

% finally, shift the points according to map if needed
if exist('map','var') && ~isempty(map)
    for ss=1:nb_grids
        S(ss).knots = map(S(ss).knots);
    end
end

% and possibly fix weights
if exist('weights_coeff','var') && ~isempty(weights_coeff)
    for ss=1:nb_grids
        S(ss).weights = S(ss).weights*weights_coeff;
    end
end

% now store the coeff value. It has to be stored after the first loop, becuase tensor_grid returns a grid
% WITHOUT coeff field, and Matlab would throw an error (Subscripted assignment between dissimilar structures)

for ss=1:nb_grids
    S(ss).coeff=coeff_condensed(ss);
end

% similarly for the multiidx generating each tensor grid
ss=1;
for j=1:nn
    if coeff(j)~=0
        i = C(j,:);
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
