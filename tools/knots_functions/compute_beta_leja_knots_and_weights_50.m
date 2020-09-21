function [X,W] = compute_beta_leja_knots_and_weights_50(alpha,beta,x_a,x_b)
%
% [X,W] = compute_beta_leja_knots_and_weights_50(alpha,beta,x_a,x_b)
% returns 50 collocation points, stored in the vector X of length 50,
% and the corresponding weights, stored in the matrix W of dim. 50x50
% (the weights for rule with P points are stored as P-th column of W),
% for the weighted Leja sequence for integration 
% w.r.t to the weight function 
%
% rho(x)=Gamma(alpha+beta+2)/ (Gamma(alpha+1)*Gamma(beta+1)*(x_b-x_a)^(alpha+beta+1)) * (x-x_a)^alpha * (x_b-x)^beta
% i.e. the density of a Beta random variable with range [x_a,x_b], alpha,beta>-1.
%
% Knots and weights are computed following the work
% 
% A. Narayan, J. Jakeman, "Adaptive Leja sparse grid constructions for stochastic collocation and high-dimensional approximation"
% SIAM Journal on Scientific Computing,  Vol. 36, No. 6, pp. A2952--A2983, 2014
%
%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

%-------------------------------------------------------
% Computation of the nodes

% Tolerance for node computation
tol = 1e-16;
% Initial and refinement step size
h0 = 1e-2;
% Initial or coarse search grid
x = -1:h0:1;
% Corresponding weight function values
w = (1-x).^(0.5*alpha) .* (1+x).^(0.5*beta);
% Initializing node vector
X = (alpha+1)/(alpha+beta+2);

% number of points to be computed
NMAX = 50; 
disp = 0:5:NMAX; % indicating if nodes function and attained minimum is plotted

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
        if ind==1 
            x_fine = x_fine(ind):h:x_fine(ind+1); % refine grid around maximum
        elseif ind == length(x_fine)
            x_fine = x_fine(ind-1):h:x_fine(ind); % refine grid around maximum
        else
            x_fine = x_fine(ind-1):h:x_fine(ind+1); % refine grid around maximum
        end
        w_fine = (1-x_fine).^(0.5*alpha) .* (1+x_fine).^(0.5*beta); % compute weights on finer grid
        % compute node function values on finer grid
        w_fine = abs(prod(repmat(x_fine,n-1,1)-repmat(X',1,length(x_fine)),1)) .* w_fine;
        % Search for maximum on finer grid
        [~,ind] = max(w_fine);    
    end

    % Update node vector 
    X(n) = x_fine(ind);

%     % uncomment this for plots of function to be minimized    
%     if (islogical(disp) && disp) || (isnumeric(disp) && ~isempty(find(disp==n)))
%        figure()
%        plot(x,w,'-',X(n),w_fine(ind),'o','LineWidth',2)
%        grid on;
%        title(sprintf('Node function and its maximum for n = %i',n))
%     end

end

% modifies points according to x_a and x_b (the weigths are unaffected)
X = ((x_a+x_b) - (x_b-x_a)*X) / 2;

%-------------------------------------------------------
% Computation of corresponding weights

tic;
W = zeros(NMAX); % we store quadrature weights for rule with P points as P-th column

for n= 1:NMAX
    nodes = X(1:n);
    [x_quad,w_quad] = knots_beta(ceil(n/2),alpha,beta,x_a,x_b);
    nnn = 1:n;
    for k=nnn
        W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
    end
end
toc;

end

