function m = lev2knots_tripling(i)

% m = lev2knots_tripling(i)
%
% relation level / number of points:
%    m = 3^{i-1}, for i>1
%    m=0          for i=0
%
% i.e. m = 1,3,9,27


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

m = 3.^(i-1);
m(i==0) = 0;

