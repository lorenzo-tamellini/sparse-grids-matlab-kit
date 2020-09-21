function [X,W] = compute_gamma_leja_knots_and_weights_50(alpha,beta)
%
% [X,W] = compute_gamma_leja_knots_and_weights_50(alpha,beta)
%
% returns 50 collocation points, stored in the vector X of length 50,
% and the corresponding weights, stored in the matrix W of dim. 50x50,
% (the weights for rule with P points are stored as P-th column of W),
% for the weighted Leja sequence for integration 
% w.r.t to the weight function 
%
% rho(x)=beta^(alpha+1)/Gamma(alpha+1)*x^alpha*exp(-beta*x) 
%
% i.e. the density of a standard Gamma random variable 
% with range [0,+inf), alpha>-1, beta>0.
%
% Knots and weights are computed following the work
% 
% A. Narayan, J. Jakeman, "Adaptive Leja sparse grid constructions for stochastic collocation and high-dimensional approximation"
% SIAM Journal on Scientific Computing,  Vol. 36, No. 6, pp. A2952--A2983, 2014
%
% The knots are computed over the interval [0,300]. 
% This choice of interval is compatible with -1<alpha<=30.
% If alpha>30 a larger interval is required and an error is raised. 
%
%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

if alpha>30
    error('SparseGKit:ToBeCoded',strcat('alpha=', num2str(alpha), ' is too large.',...
    ' Select alpha<=30'));
end

%-------------------------------------------------------
% Computation of the nodes

% Tolerance for node computation
tol = 1e-13;
% Initial and refinement step size
h0 = 1e-2;
% Initial or coarse search grid
x = 0:h0:300;
% Corresponding weight function values
w = x.^(0.5*alpha) .* exp(-0.5 * x);
% Initializing node vector
X = (alpha+1)/beta; % 0;

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
    
     while(h > tol) && (length(x_fine)>1) 
       
        h = h0*h; % refine step size
        if ind == 1
            x_fine = x_fine(ind):h:x_fine(ind+1);
        else
        x_fine = x_fine(ind-1):h:x_fine(ind+1); % refine grid around maximum
        end
        w_fine = x_fine.^(0.5*alpha) .* exp(-0.5 * x_fine); % compute weights on finer grid
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

%-------------------------------------------------------
% Computation of corresponding weights

tic;
W = zeros(NMAX); % we store quadrature weights for rule with P points as P-th column

for n = 1:NMAX
    nodes = X(1:n);
    [x_quad,w_quad] = knots_gamma(ceil(n/2),alpha,1);
    nnn = 1:n;
    for k=nnn
        W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
    end
end
toc;

% modifies points according to beta (the weigths are unaffected)
X=X/beta;

end

