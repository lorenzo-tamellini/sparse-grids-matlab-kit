function [x,w,floc]=knots_gamma_leja(n,alpha)

% [x,w] = knots_gamma_leja(n,alpha)
%
% returns the collocation points (x) and the weights (w) 
% for the weighted Leja sequence for integration 
% w.r.t to the weight function 
%
% rho(x)=x^alpha*exp(-x) 
%
% i.e. the density of a standard Gamma random variable 
% with parameters alpha>-1 and beta = 1.
%
% Knots and weights are computed (up to 50) following the work
% 
% A. Narayan, J. Jakeman, "Adaptive Leja sparse grid constructions for stochastic collocation and high-dimensional approximation"
% SIAM Journal on Scientific Computing,  Vol. 36, No. 6, pp. A2952–A2983, 2014
%
% an error is raised if more than 50 points are requested.
%
% knots are sorted increasingly before returning (weights are returned in the corresponding order)
%
% The knots are computed over the interval [x_l,x_r] = [0,450]. 
% This choice of interval is compatible with -1<=alpha<=100.
% If alpha>100 a larger interval is required and an error is raised. 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

floc = localfunctions;

if n>50
   error('SparseGKit:OutOfTable',strcat('this number of points is not available:',num2str(n)))
    
elseif alpha>100
        error('SparseGKit:OutOfTable',strcat('alpha=', num2str(alpha), ' is too large.',...
        ' Select alpha<100 or enlarge the interval [x_l, x_r] in the function get_GammaLejaKnotsAndWeights50'))

else
    
    [X,W] = get_GammaLejaKnotsAndWeights50(alpha);
    x = X(1:n);
    w = W(1:n,n);
    
    
    % sort knots increasingly and weights accordingly. Weights need to be row vectors
    [x,sorter]=sort(x);
    w=w(sorter)';

end

end

function [X,W] = get_GammaLejaKnotsAndWeights50(alpha)

    %-------------------------------------------------------
    % Computation of the nodes

    % Tolerance for node computation
    tol = 1e-16;
    % Initial and refinement step size
    h0 = 1e-2;
    % Initial or coarse search grid
    x_l = 0; x_r = 450;
    x = x_l:h0:x_r;
    % Corresponding weight function values
    w = x.^(0.5*alpha) .* exp(-0.5 * x);
    % Initializing node vector
    X = 0;

    % number of points to be computed
    NMAX = 50; 
    % disp = 0:5:NMAX; % indicating if nodes function and attained minimum is plotted

    % Loop over length of node vector
    for n = 2:NMAX
        % Update weighted node function to be maximized
        x_k = X(n-1);
        w = abs(x-x_k) .* w;

        % Search for maximum on coarse level
        [~,ind] = max(w);

        % Adaptively refine search grid around maximum on current grid
        h = h0;
        x_fine = x; 
        while(h > tol)
            h = h0*h; % refine step size
            x_fine = x_fine(ind-1):h:x_fine(ind+1); % refine grid around maximum
            w_fine = x_fine.^(0.5*alpha) .* exp(-0.5 * x_fine); % compute weights on finer grid
            % compute node function values on finer grid
            w_fine = abs(prod(repmat(x_fine,n-1,1)-repmat(X',1,length(x_fine)),1)) .* w_fine;
            % Search for maximum on finer grid
            [~,ind] = max(w_fine);    
        end

        % Update node vector 
        X(n) = x_fine(ind);

    end
    
    %-------------------------------------------------------
    % Computation of corresponding weights

    tic;
    W = zeros(NMAX); % we store quadrature weights for rule with P points as P-th column

    for n= 1:NMAX
        nodes = X(1:n);
        [x_quad,w_quad] = knots_gamma(ceil(n/2),alpha,1);
        nnn = 1:n;
        for k=nnn
            W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
        end
    end
    toc;

end

