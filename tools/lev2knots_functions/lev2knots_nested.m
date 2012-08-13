% ========================================================
%   relation level / number of points:
%    m = 2^{i-1}+1, for i>1
%    m=1            for i=1
%    m=0            for i=0
%   when using Clenshaw-Curtis or Gauss-Patterson formulae,
%   gives nested sequences of points
% ========================================================
%
%   [m] = lev2knots_nested(i)
%   i: level in each direction
%   m: number of points to be used in each direction

function [m] = lev2knots_nested(i)

m = 2.^(i-1)+1;
for k=1:length(m(:))
    if i(k)==1 
        m(k) = 1;
    end
    if i(k)==0
        m(k) = 0;
    end
end

