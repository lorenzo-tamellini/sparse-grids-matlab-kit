function [x,w] = knots_beta_leja(n,alpha,beta,x_a,x_b,type,saving_flag)
%
% [x,w] = KNOTS_BETA_LEJA(n,alpha,beta,x_a,x_b,<type>,<saving-flag>)
%
% returns the first n collocation points (x) and the corresponding weights (w), 
% for the weighted Leja sequence for integration w.r.t to the weight function 
%
% rho(x)=Gamma(alpha+beta+2)/ (Gamma(alpha+1)*Gamma(beta+1)*(x_b-x_a)^(alpha+beta+1)) * (x-x_a)^alpha * (x_b-x)^beta
% i.e. the density of a Beta random variable with range [x_a,x_b], alpha,beta>-1.
%
% Knots and weights are computed following the work
% 
% A. Narayan, J. Jakeman, "Adaptive Leja sparse grid constructions for stochastic collocation and high-dimensional approximation"
% SIAM Journal on Scientific Computing,  Vol. 36, No. 6, pp. A2952--A2983, 2014
%
% Knots are sorted increasingly before returning (weights are returned in the corresponding order).
% 
% if <saving_flag> is set to 'on_file', the function will compute max(n,50) knots and weights, store them in a .mat file
% in the local folder and load results at next calls with the same values of alpha and beta. Note that nodes on
% different intervals [x_a x_b] can be obtained by linear translations from the reference interval [-1 1],
% therefore calls like
%
% knots_beta_leja(n,alpha,beta,a1,b1,<type>,'on_file')
% knots_beta_leja(n,alpha,beta,a2,b2,<type>,'on_file')
%
% can use the same precomputed file. Load command will look for the file in all folders added to the path of matlab
%
% if no saving flag is provided, i.e., the function is called as knots_beta_leja(n,alpha,beta,x_a,x_b,<type>), 
% computations of nodes and weights will not be saved. This might become time-consuming for large sparse grids
% (the computation of the Leja knots is done by brute-force optimization).
% 
% Follows the description of the choices of <type>.
%
% -----------------------------------------------------------
%
% [X,W] = KNOTS_BETA_LEJA(n,alpha,beta,x_a,x_b,'line') given X(1)=(alpha+1)/(alpha+beta+2)
%  recursively defines the n-th point by 
%
%   X_n= argmax_[x_a,x_b] prod_{k=1}^{n-1} abs(X-X_k)
% 
%
% [X,W] = KNOTS_BETA_LEJA(N,mi,sigma,'sym_line') given X(1)=0 recursively
%   defines the n-th and (n+1)-th point by 
%
%   X_n= argmax prod_{k=1}^{n-1} abs(X-X_k)
%   X_(n+1) = symmetric point of X_n with respect to 0

%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------


if nargin == 7 
    switch saving_flag
        case 'on_file'
            save_and_load = true;
            filename = ['SGK_knots_weights_beta_',type,'_leja_file_alpha_',num2str(alpha),'_beta_',num2str(beta),'.mat'];
        otherwise
            error('SparseGKit:WrongInput','unknown saving option for knots_beta_leja')
    end
else
    save_and_load = false;
end

switch type 
    %--------------------------------------------------
    case 'line'
        if save_and_load
            if exist(filename,'file') == 2 
               % load from file
               load(filename,'Xall','Wall'); 
               if n > length(Xall)
                   error('SparseGKit:OutOfTable',['not enough beta points computed on file. Remove file ',filename,...
                       ' and run again KNOTS_BETA_LEJA with ''on_file'' option and larger n'])
               end
               X = Xall(1:n);
               W = Wall(1:n,n);
            else
               % compute a large number of nodes
               [Xall,Wall] = get_BetaLejaKnotsAndWeights(max(50,n),alpha,beta,'all');
               % save on file     
               save(filename,'Xall','Wall');  
               X = Xall(1:n);
               W = Wall(1:n,n);
            end
        else
            % just compute what asked
            [X,W] = get_BetaLejaKnotsAndWeights(n,alpha,beta,'current');
        end
        
    %--------------------------------------------------    
    case 'sym_line'
        if alpha ~= beta
            warning('SparseGKit:InconsistentInput','The shape parameters alpha and beta are not equal. Hence, the beta pdf is not symmetric and working with symmetric knots is not recommended.'); 
        end
        if save_and_load
            if exist(filename,'file') == 2 
                % load from file
                load(filename,'Xall','Wall');
                if n > length(Xall)
                    error('SparseGKit:OutOfTable',['not enough beta points computed on file. Remove files ',filename,...
                        ' and run again KNOTS_BETA_LEJA with ''on_file'' option and larger n'])
                end
                X = Xall(1:n);
                W = Wall(1:n,n);                
            else
                % disp('compute 50 nodes')
                [Xall,Wall] = get_SymBetaLejaKnotsAndWeights(max(50,n),alpha,beta,'all');
                % save on file
                save(filename,'Xall','Wall');
                X = Xall(1:n);
                W = Wall(1:n,n);
            end
        else
            % just compute what asked
            [X,W] = get_SymBetaLejaKnotsAndWeights(n,alpha,beta,'current');
        end        
        
    %--------------------------------------------------    
    otherwise
        error('SparseGKit:WrongInput','unknown Leja type')
end


% modifies points according to x_a and x_b 
X = ((x_a+x_b) + (x_b-x_a)*X) / 2; 

% sort knots increasingly and weights accordingly. Weights need to be row vectors
[x,sorter]=sort(X);
w=W(sorter)';

end



function [X,W] = get_BetaLejaKnotsAndWeights(n,alpha,beta,which_weights) % x_a=-1,x_b=1
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
            [x_quad,w_quad] = knots_beta(ceil(n/2),alpha,beta,-1,1);
            nnn = 1:n;
            for k=nnn
                W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
            end
        end
        
    case 'current' % only weights for formula with NMAX nodes
        W = zeros(NMAX,1);
        nodes = X;
        [x_quad,w_quad] = knots_beta(ceil(NMAX/2),alpha,beta,-1,1);
        nnn = 1:NMAX;
        for k=nnn
            W(k) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
        end
end

end



function [X,W] = get_SymBetaLejaKnotsAndWeights(n,alpha,beta,which_weights)
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

% number of points to be computed
NMAX = n; 

% Initializing node vector
X = 0; % first node 

% Loop over length of node vector
for j = 2:NMAX
    % Update weighted node function to be maximized
    x_k = X(j-1);
    w = abs(x-x_k) .* w;

    if mod(j,2)==1
        % Symmetric node
        X(j) = -X(j-1);
        
    else

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
            w_fine = abs(prod(repmat(x_fine,j-1,1)-repmat(X',1,length(x_fine)),1)) .* w_fine;
            % Search for maximum on finer grid
            [~,ind] = max(w_fine);    
        end

        % Update node vector 
        X(j) = x_fine(ind);

    end
end

%-------------------------------------------------------
% Computation of corresponding weights

switch which_weights
    
    case 'all' % all weights for formulas with 1 node up to NMAX nodes
        W = zeros(NMAX); % we store quadrature weights for rule with P points as P-th column        
        for n = 1:NMAX
            nodes = X(1:n);
            [x_quad,w_quad] = knots_beta(ceil(n/2),alpha,beta,-1,1);
            nnn = 1:n;
            for k=nnn
                W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
            end
        end      
        
    case 'current'  % only weights for formula with NMAX nodes
        W = zeros(NMAX,1);
        nodes = X;
        [x_quad,w_quad] = knots_beta(ceil(NMAX/2),alpha,beta,-1,1);
        nnn = 1:NMAX;
        for k=nnn
            W(k) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
        end
end
end

