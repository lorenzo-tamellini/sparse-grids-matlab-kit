function [S,C] = smolyak_grid_multiindeces(C,knots,lev2knots,map,weights_coeff)

% SMOLYAK_GRID_MULTIINDICES produces a sparse grid starting from a multiindex-set rather than
% from a rule IDXSET(I) <= W.
% 
% [S,C] = SMOLYAK_GRID_MULTIINDICES(C,KNOTS,LEV2KNOTS) uses the multiindex set C. C must be
%       in lexicographic order and admissible. 
%
% [S,C] = SMOLYAK_GRID_MULTIINDICES(C,KNOTS,LEV2KNOTS,MAP,WEIGHTS_COEFF) can be used as an alternative
%       to generate a sparse grid on a hyper-rectangle, 
%
% See also CHECK_SET_ADMISSIBILITY for admissibility check, and SMOLYAK_GRID for further information on 
% KNOTS, LEV2KNOTS, MAP, WEIGHTS_COEFF and on the sparse grid data structure S


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
% naif implementation; exploit partial ordering of the sequence of
% multiindices

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

% disp(strcat('using:',num2str(nn),'multiindices'))


coeff=ones(1,nn); % initialize coefficients to 1: all c survive

for i=1:nn % scroll c
    for j=i+1:nn % scroll c2, the following rows
        d=C(j,:)-C(i,:);
        if d<=1 & d>=0
            coeff(i)=coeff(i) + (-1)^sum(d);
        end
    end

end

% now for those c who survived, compute the grid, and multiply with appropriate coefficients

for j=1:nn
    if coeff(j)~=0
        i = C(j,:);       % level in each direction
        m = apply_lev2knots(i,lev2knots,N); % n. of points in each direction
        S(j) = tensor_grid(N,m,knots);
        S(j).weights=S(j).weights*coeff(j);
    end
end




% finally, shift the points according to map if needed
if exist('map','var')
    for j=1:nn
        if coeff(j)~=0
            S(j).knots = map(S(j).knots);
        end
    end
end

% and possibly fix weights
if exist('weights_coeff','var')
    for j=1:nn
        if coeff(j)~=0
            S(j).weights = S(j).weights*weights_coeff;
        end
    end
end



% now store the coeff value. It has to be stored after the first loop, becuase tensor_grid returns a grid
% WITHOUT coeff field, so that if you add it then you get a mismatch between the fields in output and the fields
% in the struct array you are using.

for j=1:nn
    if coeff(j)~=0
        S(j).coeff=coeff(j);
    end
end


end





function m = apply_lev2knots(i,lev2knots,N)
    
% N could be deduced by N but it's better passed as an input, to speed computation
% init m to zero vector
m=0*i;

% next, iterate on each direction 1,...,N. 

for n=1:N
    m(n) = lev2knots{n}(i(n));
end

end