function isl = islexico_tol(a,b,TOL)


% ISLEXICO_TOL returns true if two vectors A, B and returns true if A<B lexicogrphically up to numerical noise of size TOL:
%
% ii = [1 3];
% jj = [1 4]; % islexico(ii,jj) -> true
% 
% ii = [1+1e-13 3];
% islexico_tol(ii,jj,1e-12) -> true
% islexico_tol(ii,jj,1e-14) -> false

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


[~,idx]=find(abs(a-b)>TOL,1); % find the pos of the first element of |a-b| larger than TOL. if none, idx is empty

% return true if v is empty or if the first component for which a and b differ more than TOL is smaller in a
% than in b
isl = isempty(idx) || a(idx) - b(idx) <0 ;

end