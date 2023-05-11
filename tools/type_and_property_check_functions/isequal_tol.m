function iseq = isequal_tol(ii,jj,TOL)


% ISEQUAL_TOL compares two vectors and returns true if they are component-wise less than TOL apart:
%
% ii    = [3 5 1];
% noise = 1e-12*randn(size(ii));
% jj = ii + noise;
%
% isequal_tol(ii,jj,1e-11) -> TRUE
% isequal_tol(ii,jj,1e-14) -> FALSE


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


iseq = max(abs(ii-jj))<TOL;
    
end
