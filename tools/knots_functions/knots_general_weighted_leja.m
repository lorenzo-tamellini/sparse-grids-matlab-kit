function [x,w] = knots_general_weighted_leja(n,X,W)

% [x,w] = knots_general_weighted_leja(n,X,W)
%
% returns the first n collocation points (x) and weights (w), 
% from the vector X, collecting the first 50 weighted Leja points, 
% and the matrix W of dim. 50x50, collecting the corresponing weights
% (the weights for rule with P points are stored as P-th column of W), respectively.  
% 
% The inputs X and W are generated either by 
% compute_gamma_leja_knots_and_weights_50 or compute_beta_leja_knots_and_weights_50
% for the integration w.r.t. Gamma or Beta densities, respectively. 
%
% An error is raised if more than 50 points are requested.
%
% Knots are sorted increasingly before returning (weights are returned in
% the corresponding order).
%
%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

if n>50
   error('SparseGKit:OutOfTable',strcat('this number of points is not available:',num2str(n)))
    
else
    
    x = X(1:n);
    w = W(1:n,n);
        
    % sort knots increasingly and weights accordingly. Weights need to be row vectors
    [x,sorter]=sort(x);
    w=w(sorter)';

end

end

