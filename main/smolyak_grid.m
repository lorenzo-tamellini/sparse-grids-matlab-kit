% ============================================================
%  Generates a Smolyak sparse grid (and corresponding quadrature weights)
%  as a linear combination of full tensor grids, employing formula (2.9)
%  of [Nobile-Tempone-Webster, SINUM 46/5, pages 2309-2345] 
%
%  The sparse grid information is stored as a vector of "tensor grids", 
%  each "tensor grid" S(j) is a three field structure:
%    S(j).knots: vector containing the tensor grid knots
%    S(j).weights: vector containing the corresponding weights
%    S(j).size: size of the tensor grid = prod(m)
%    S(j).coeff: how many times the tensor grid appears in the sparse grid (with sign)
%  the index j runs over all points in the level set Y(w,N) 
%
%  usage:
%   [S,C] = smolyak_grid(N,w,knots,lev2knots,idxset)
%   input
%     N: dimension
%     w: level in the Smolyak formula (w>=0)
%     knots: function defining the 1D gauss knots
%        header: [x,w]=knots(m)  (generates m knots and weights)
%     lev2knots: function defining the relation between level and 
%        numebr of knots.   header: m=lev2knots(i)
%     idxset (optional): function defining the index set for the levels
%        a multiindex is in the set if idxset(i)<=w (default is sum(i-1))
%   output
%     S: structure containing the information on the sparse grid 
%        (vector of tensor grids; see above)
%     C: multi-index set used to generate the sparse grid  
% ============================================================

function [S,C] = smolyak_grid(N,w,knots,lev2knots,idxset)

if nargin==4
    idxset=@(i) sum(i-1);
end

if w==0
    i = ones(1,N);
    m = lev2knots(i);
    S(1) = tensor_grid(N,m,knots);
    S(1).coeff=1;
    disp('using 1 multiindex')
    C=i;
else
    
    % build the list of multiindices in the set: idxset(i)<=w
    C=multiidx_gen(N,idxset,w,1);
    % disp(['using ',num2str(size(C,1)),' multiindeces'])

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
    % C is partially ordered (lexicographis): [x x x x i], are listed increasing with i, 
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
            m = lev2knots(i); % n. of points in each direction
            S(j) = tensor_grid(N,m,knots);
            S(j).weights=S(j).weights*coeff(j);
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
