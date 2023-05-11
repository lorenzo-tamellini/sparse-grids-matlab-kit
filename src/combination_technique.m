function coeff = combination_technique(I)

% COMBINATION TECHNIQUE computes the coefficients of the combination technique expression of a sparse grid
%
% COEFF = COMBINATION_TECHNIQUE(I) takes as input a multiidx set I, each associated to a Delta operator
%       in the sparse grids construction. It then computes the coefficients of the combination technique
%       expression of the same sparse grids, i.e.,  the same sparse grids computed as linear combination
%       of tensor grids. 
%
%       I must be a multiidx set matrix (one row per index), admissible and lexicographically sorted.
%       These two properties will **not** be checked by the function.
%
%       COEFF is a vector containing the combination technique coefficients

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

nn=size(I,1);
coeff=ones(1,nn); % initialize coefficients to 1: all c survive



% I can at least restrict the search to multiindices whose first component is c(i) + 2, so I define
[~,bookmarks]=unique(I(:,1),'first');
bk = [bookmarks(3:end)'-1 nn nn];
% i.e. those who begin with 1 end at bookmark(3)-1, those who begin with 2-1 end at bookmark(4) and so on,
% until there's no multiindex with c(i)+2
%
% an old piece of code I still want to have written somewhere
% range = i+ find( C(i+1:end,1)==cc(1)+2, 1 ) - 1;


for i=1:nn % scroll c
    cc = I(i,:);
    % recover the range in which we have to look. Observe that the first column of C contains necessarily 1,2,3 ...
    % so we can use them to access bk
    range=bk(cc(1));
    for j=i+1:range
        % scroll c2, the following rows
        d=I(j,:)-cc;
        if max(d)<=1 && min(d)>=0  % much faster to check then if only 0 and 1 appears. Also, && is short-circuited,
            % so if max(d)<=1 is false the other condition is not even checked
            coeff(i)=coeff(i) + (-1)^sum(d);
        end
    end
end
