function [x,w] = knots_gamma_leja(n,alpha,beta,saving_flag)
%
% [x,w] = knots_gamma_leja(n,alpha,beta,<saving-flag>)
%
% returns the first n collocation points (x) and the corresponding weights (w), 
% for the weighted Leja sequence for integration 
% w.r.t the weight function 
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
% The maximum number of points that can be computed is 50, for compatibility with the interval [0,300].
%
% if <saving_flag> is set to 'on_file', the function will compute 50 knots and weights, store them in a .mat file
% in the local folder and load results at next calls with the same values of alpha. Note that different values of beta
% act only as rescaling factor, therefore the calls
%
% knots_gamma_leja(n,alpha,beta1,'on_file')
% knots_gamma_leja(n,alpha,beta2,'on_file')
%
% can use the same precomputed file. Load command will look for the file in all folders added to the path of matlab
%
% if no saving flag is provided, i.e., the function is called as knots_gamma_leja(n,alpha,beta), 
% computations of nodes and weights will not be saved. This might become time-consuming for large sparse grids
% (the computation of the Leja knots is done by brute-force optimization).
%
%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

if alpha>30
    error('SparseGKit:ToBeCoded',strcat('alpha=', num2str(alpha), ' is too large. Select alpha<=30'));
end

if n>50
   error('SparseGKit:OutOfTable',strcat('this number of points is not available:',num2str(n)))
end


if nargin == 4 
    switch saving_flag
        case 'on_file'
            save_and_load = true;
            filename = ['SGK_knots_weights_gamma_leja_file_alpha_',num2str(alpha),'.mat'];
        otherwise
            error('SparseGKit:WrongInput','unknown saving option for knots_gamma_leja')
    end
else
    save_and_load = false;
end



if save_and_load
    if exist(filename,'file') == 2
        % load from file
        load(filename,'Xall','Wall');
        if n > length(Xall) % this is probably useless since we hard-coded a limit on n=50 and the next case of this if computes 50 nodes ... 
                            % but is good to keep it in case we relax the n=50 limit
            error('SparseGKit:OutOfTable',['not enough gamma points computed on file. Remove files ',filename_knots,' and ',filename_weights,...
                ' and run again KNOTS_GAMMA_LEJA with ''on_file'' option and larger n'])
        end
        X = Xall(1:n);
        W = Wall(1:n,n);
    else
        % compute 50 nodes
        [Xall,Wall] = get_GammaLejaKnotsAndWeights(50,alpha,'all'); % here beta=1
        save(filename,'Xall','Wall');
        X = Xall(1:n);
        W = Wall(1:n,n);
    end
else
    % just compute
    [X,W] = get_GammaLejaKnotsAndWeights(n,alpha,'current'); % here beta=1
end


% modifies points according to beta 
X=X/beta;

% sort knots increasingly and weights accordingly. Weights need to be row vectors
[x,sorter]=sort(X);
w=W(sorter)';

end

function [X,W] = get_GammaLejaKnotsAndWeights(n,alpha,which_weights) % here beta fixed = 1
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
X = (alpha+1); 

% number of points to be computed
NMAX = n; 

% Loop over length of node vector
for j = 2:NMAX
    
    % Update weighted node function to be maximized
    x_k = X(j-1);
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
        w_fine = abs(prod(repmat(x_fine,j-1,1)-repmat(X',1,length(x_fine)),1)) .* w_fine;
        
        % Search for maximum on finer grid
        [~,ind] = max(w_fine); 

    end

    % Update node vector 
    X(j) = x_fine(ind);
    
end   

    
%-------------------------------------------------------
% Computation of corresponding weights
switch which_weights
    
    case 'all' % all weights for formulas with 1 node up to NMAX nodes
        W = zeros(NMAX); % we store quadrature weights for rule with P points as P-th column        
        for n = 1:NMAX
            nodes = X(1:n);
            [x_quad,w_quad] = knots_gamma(ceil(n/2),alpha,1); % here 1 = beta
            nnn = 1:n;
            for k=nnn
                W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
            end
        end
        
    case 'current' % only weights for formula with NMAX nodes
        W = zeros(NMAX,1);
        nodes = X;
        [x_quad,w_quad] = knots_gamma(ceil(NMAX/2),alpha,1); % here 1 = beta
        nnn = 1:NMAX;
        for k=nnn
            W(k) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
        end
end
 
end
