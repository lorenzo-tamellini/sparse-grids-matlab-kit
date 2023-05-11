%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, B. Sprungk
% See LICENSE.txt for license
%----------------------------------------------------



%% Computation of Normal Leja Nodes and weights
% for the weighted Leja sequence for integration 
% w.r.t to the standard gaussian weight function 
% rho(x)=1/sqrt(2*pi)*exp(-x^2/2) 

clear

tic;
% Tolerance for node computation
tol = 1e-16;
% Initial and refinement step size
h0 = 1e-2;
% Initial or coarse search grid
x = -40:h0:40;
% Corresponding weight function values
w = exp(-0.25 * x.^2);
% Initializing node vector
X = 0;

% number of points to be computed
NMAX = 150; 
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
        x_fine = x_fine(ind-1):h:x_fine(ind+1); % refine grid around maximum
        w_fine = exp(-0.25 * x_fine.^2); % compute weights on finer grid
        % compute node function values on finer grid
        w_fine = abs(prod(repmat(x_fine,n-1,1)-repmat(X',1,length(x_fine)),1)) .* w_fine;
        % Search for maximum on finer grid
        [~,ind] = max(w_fine);    
    end
    
    % Update node vector 
    X(n) = x_fine(ind);

%   uncomment this for plots of function to be minimized    
%     if (islogical(disp) && disp) || (isnumeric(disp) && ~isempty(find(disp==n))),
%        figure()
%        plot(x,w,'-',X(n),w_fine(ind),'o','LineWidth',2)
%        grid on;
%        title(sprintf('Node function and its maximum for n = %i',n))
%     end
end
toc;


%-------------------------------------------------------
% Computation of corresponding weights

tic;
W = zeros(NMAX); % we store quadrature weights for rule with P points as P-th column


for n= 1:NMAX
    nodes = X(1:n);
    [x_quad,w_quad] = knots_normal(ceil(n/2),0,1);
    nnn = 1:n;
    for k=nnn
        W(k,n) = dot(lagr_eval(nodes(k), nodes(nnn~=k),x_quad), w_quad);
    end
end
toc;

% plot(X) when X is matrix plots the columns of X, therefore

figure
plot(W,'-') % plots weights of each quadrature rule
semilogy(abs(W),'-') % plots weights of each quadrature rule

figure
plot(W','-') % plots trend of weights of each quadrature rule

% save('GaussianLejaPrecomputedKnotsAndWeights','X','W')
% save('gaussian_leja.txt','X','W','-ascii','-double')


%% Testing each Normal-Leja quadrature on computation of moments of a gaussian random variable

clear

% Moments of standard normal distribution
p_max = 51;
mom = zeros(1,1+p_max);
for p = 0:2:p_max
    mom(p+1) = prod(1:2:p);
end

% Quadrature error for polynomials
imax = 50;
rel_err = zeros(1+p_max,imax);
rel_errGH = zeros(1+p_max,imax);

% for each formula, we test its approximation of increasing moments

disp('plotting relative error for several quadrature rules just generated - press ENTER to move to the next plot or CTRL+C to stop')
for n=10:5:50
    % Normal-Leja quadrature rule using n nodes
    [x_Lj,w_Lj] = knots_normal_leja(n,0,1,'line');
    
    % Gauss-Hermite quadrature of same accuracy
    [x_GH,w_GH] = knots_normal(ceil(n/2),0,1);
    
    for p=0:p_max
        rel_err(1+p,n) = (abs(mom(1+p) - dot(x_Lj.^p,w_Lj) ))./mom(1+p);
        rel_errGH(1+p,n) = (abs(mom(1+p) - dot(x_GH.^p,w_GH) ))./mom(1+p);
    end
    
    % Plotting errors
    figure
    semilogy(0:p_max, rel_err(:,n),'o',0:p_max, rel_errGH(:,n),'+','LineWidth',2)
    grid on
    ylim([1e-20 1e5])
    set(gca,'FontSize',14)
    legend(sprintf('Normal-Leja (n=%i)',n),sprintf('Gauss-Hermite (n=%i)',ceil(n/2)))
    set(legend,'Location','SouthOutside')
    xlabel('Polynomial degree p')
    ylabel('Relative quadrature error')
    title('paused - press ENTER to move to the next plot or CTRL+C to stop','FontSize',10)
    
    pause
    
    close
end



%% now we do the opposite: fix function and increase point in quad rule. 
% The Leja formula is built to be optimal in interpolation, so it does not perform
% very well in quadrature, as mentioned in the reference paper by A. Narayan and J. Jakeman

clear

imax=38;

% function to be integrated 

%f=@(x) 1./(1+0.1*x.^2); 
f=@(x) 1./(2+exp(0.2*x)); 
%f=@(x) cos(3*x+1); 

quad_Lj=zeros(1,imax);
quad_GH=zeros(1,imax);
nb_pts =zeros(1,imax);

% refining quad rule
for i=1:imax
    
    n = lev2knots_lin(i);
    nb_pts(i) = n;

    [x_Lj,w_Lj]=knots_normal_leja(n,0,1,'line');   
    [x_GH,w_GH] = knots_normal(n,0,1);
    
    quad_Lj(i) = dot(f(x_Lj),w_Lj);
    quad_GH(i) = dot(f(x_GH),w_GH);
end

% repeat for GK knots, that have a different limit of precomputation
quad_GK=zeros(1,5);
nb_pts_GK =zeros(1,5);

for i=1:5
    
    n = lev2knots_GK(i);
    nb_pts_GK(i) = n;

    [x_GK,w_GK] = knots_GK(n,0,1);
    
    quad_GK(i) = dot(f(x_GK),w_GK);
    
end

% exact integral
[x_GH,w_GH] = knots_normal(160,0,1);
exact = dot(f(x_GH),w_GH);
err_Lj=abs(quad_Lj - exact);
err_GH=abs(quad_GH - exact);
err_GK=abs(quad_GK - exact);


% Plotting errors
figure
semilogy(nb_pts, err_Lj,'-xr','LineWidth',2,'DisplayName','Normal Leja pts')
grid on
hold on
semilogy(nb_pts, err_GH,'-ob','LineWidth',2,'DisplayName','Gauss Hermite pts')
semilogy(nb_pts_GK, err_GK,'-^k','LineWidth',2,'DisplayName','GK pts')


legend show
set(legend,'Location','SouthWest')
%set(gca,'FontSize',14)

ylim([1e-16 10])


%% because in 1D Weighted Leja suffer against GK and GH, it is interesting to check 
% the effect in the multi-variate case. Here the "granularity" of Leja should pay off.
% However, the fact that Leja are not designed for quadrature still shows up

clear

% dimension of space
N=2;

% we use a simple TD rule, up to this level
w_max=15;

% function to be integrated
%f=@(x) 1./(1+0.1*norm(x).^2); 
f=@(x) 1./(2+exp(0.2*sum(x))); 
%f=@(x) 1./(2+exp(sum(x))); 

knots_Lj = @(n) knots_normal_leja(n,0,1,'line');   
knots_GH = @(n) knots_normal(n,0,1);

quad_Lj=zeros(1,w_max);
quad_GH=zeros(1,w_max);

nb_pts_Lj =zeros(1,w_max);
nb_pts_GH =zeros(1,w_max);

% we introduce auxiliary containers to recycle previous evaluations and speed up the computation
S_Lj_old=[];
Sr_Lj_old=[];
evals_Lj_old=[];

S_GH_old=[];
Sr_GH_old=[];
evals_GH_old=[];

% the convergence loop for Leja and Gauss Hermite
for w=1:w_max       
    
    disp(w); %#ok<*NOEFF>
    
    disp('Leja');
    S_Lj = create_sparse_grid(N,w,knots_Lj,@lev2knots_2step, @(i) sum(i-1), S_Lj_old); %using 2step rule to ramp faster to rules with high number of points
    Sr_Lj = reduce_sparse_grid(S_Lj);
    [res,evals]= quadrature_on_sparse_grid(f, S_Lj, Sr_Lj, evals_Lj_old, S_Lj_old, Sr_Lj_old);
    quad_Lj(w) = res;
    evals_Lj_old = evals;
    S_Lj_old=S_Lj;
    Sr_Lj_old = Sr_Lj;
    nb_pts_Lj(w) = Sr_Lj.size;

    disp('Gauss Hermite');
    S_GH = create_sparse_grid(N,w,knots_GH,@lev2knots_lin, @(i) sum(i-1), S_GH_old);
    Sr_GH = reduce_sparse_grid(S_GH);    
    [res, evals]  = quadrature_on_sparse_grid(f, S_GH, Sr_GH, evals_GH_old, S_GH_old, Sr_GH_old);
    quad_GH(w) = res;
    evals_GH_old = evals;
    S_GH_old = S_GH;
    Sr_GH_old = Sr_GH;    
    nb_pts_GH(w) = Sr_GH.size;
end


% repeat for GK 
%knots_GK = @(n) knots_GK(n);   

w_max_GK=4;

knots_GK = @(n) knots_GK(n,0,1);   

quad_GK=zeros(1,w_max_GK);
nb_pts_GK =zeros(1,w_max_GK);

S_GK_old=[];
Sr_GK_old=[];
evals_GK_old=[];


for w=1:w_max_GK       
    
    disp(w);
    disp('GK knots');
    S_GK = create_sparse_grid(N,w,knots_GK,@lev2knots_GK, @(i) sum(i-1), S_GK_old);
    Sr_GK = reduce_sparse_grid(S_GK);
    [res, evals] = quadrature_on_sparse_grid(f,S_GK,Sr_GK,evals_GK_old,S_GK_old,Sr_GK_old);
    quad_GK(w) = res;
    evals_GK_old = evals;
    S_GK_old = S_GK;
    Sr_GK_old = Sr_GK;    
    nb_pts_GK(w) = Sr_GK.size;
end


% exact integral
disp('computing reference solution');
S_GH = create_sparse_grid(N,w_max+2,knots_GH,@lev2knots_lin, @(i) sum(i-1), S_GH_old);
Sr_GH = reduce_sparse_grid(S_GH);
exact = quadrature_on_sparse_grid(f,S_GH,Sr_GH,evals_GH_old,S_GH_old,Sr_GH_old);

% errors
err_Lj=abs(quad_Lj - exact);%./abs(quad_GH(end));
err_GH=abs(quad_GH - exact);%./abs(quad_GH(end));
err_GK=abs(quad_GK - exact);%./abs(quad_GH(end));


figure
loglog(nb_pts_Lj, err_Lj,'-xr','LineWidth',2,'DisplayName','Normal Leja pts')
grid on
hold on
loglog(nb_pts_GH, err_GH,'-ob','LineWidth',2,'DisplayName','Gauss Hermite pts')
loglog(nb_pts_GK, err_GK,'-^k','LineWidth',2,'DisplayName','GK pts')


legend show
set(legend,'Location','SouthWest')



%% now we move to interpolation error in 1D. Here the fact that Leja points are constructed for interpolation will pay

clear


% function to be interpoled
f=@(x) 1./(1+0.1*x.^2); 
%f=@(x) 1./(2+exp(x)); 
%f=@(x) cos(3*x+1); 

% we sample the function on a random sample 
samplesize = 100;
sampleset = randn(1,samplesize);
F_samples = f(sampleset);

% the convergence analysis
imax=150;

err_Lj=zeros(1,imax);
err_GH=zeros(1,imax);

nb_pts =zeros(1,imax);

for i=1:imax
    
    n = lev2knots_lin(i);
    nb_pts(i) = n;
    nnn = 1:n;
        
    % here we build the lagrange interpolant for Leja and evaluate the error
    [x_Lj,w_Lj]=knots_normal_leja(n,0,1,'line');       
    interp_Lj = zeros(1,samplesize);
    for k=nnn
        interp_Lj =  interp_Lj + f(x_Lj(k))*lagr_eval(x_Lj(k), x_Lj(nnn~=k),sampleset);
    end    
    err_Lj(i) = max(abs(F_samples - interp_Lj)); 


    % repeat for Gauss Hermite
    [x_GH,w_GH]=knots_normal(n,0,1);
    interp_GH = zeros(1,samplesize);
    for k=nnn
        interp_GH =  interp_GH + f(x_GH(k))*lagr_eval(x_GH(k), x_GH(nnn~=k),sampleset);
    end    
    err_GH(i) = max(abs(F_samples - interp_GH)); 


end

% repeat for GK
imax_GK=5;
err_GK=zeros(1,imax_GK);
nb_pts_GK =zeros(1,imax_GK);

for i=1:imax_GK
    
    n = lev2knots_GK(i);
    nb_pts_GK(i) = n;
    nnn = 1:n;
    
    [x_GK,w_GK]=knots_GK(n,0,1);       
    interp_GK = zeros(1,samplesize);
    for k=nnn
        interp_GK =  interp_GK + f(x_GK(k))*lagr_eval(x_GK(k), x_GK(nnn~=k),sampleset);
    end    
    err_GK(i) = max(abs(F_samples - interp_GK)); 

end



% Plotting errors
figure
semilogy(nb_pts, err_Lj,'-xr','LineWidth',2,'DisplayName','Normal Leja pts')
hold on
semilogy(nb_pts, err_GH,'-ob','LineWidth',2,'DisplayName','Gauss Hermite pts')
semilogy(nb_pts_GK, err_GK,'-^k','LineWidth',2,'DisplayName','GK pts')
grid on

legend show
set(legend,'Location','SouthWest')


%% repeat in the multivariate case. 


clear

% the function to be interpolated 

N=2;

%f=@(x) 1./(1+0.1*sum(x.^2)); 
f=@(x) 1./(2+exp(0.2*sum(x))); 

% we evaluate the function over a sample set
samplesize = 1000;
sampleset = randn(N,samplesize);
f_sampled = f(sampleset);

% as before, we use simple TD sparse grids, up to this level
w_max=15;


% the convergence loop, with auxiliary containers to recycle evaluations between iterations
knots_Lj = @(n) knots_normal_leja(n,0,1,'line');   
knots_GH = @(n) knots_normal(n,0,1);

err_Lj=zeros(1,w_max);
err_GH=zeros(1,w_max);

nb_pts_Lj =zeros(1,w_max);
nb_pts_GH =zeros(1,w_max);

S_Lj_old=[];
Sr_Lj_old=[];
evals_Lj_old=[];

S_GH_old=[];
Sr_GH_old=[];
evals_GH_old=[];


for w=1:w_max       
    
    disp(w); %#ok<*NOEFF>
    
    disp('Leja');
    S_Lj = create_sparse_grid(N,w,knots_Lj,@lev2knots_2step, @(i) sum(i-1), S_Lj_old); %using 2step rule to ramp faster to rules with high number of points
    Sr_Lj = reduce_sparse_grid(S_Lj);
    evals_Lj = evaluate_on_sparse_grid(f, S_Lj, Sr_Lj, evals_Lj_old, S_Lj_old, Sr_Lj_old);
    err_Lj(w) = max(abs(f_sampled - interpolate_on_sparse_grid(S_Lj,Sr_Lj,evals_Lj,sampleset)));
    evals_Lj_old = evals_Lj;
    S_Lj_old=S_Lj;
    Sr_Lj_old = Sr_Lj;
    nb_pts_Lj(w) = Sr_Lj.size;

    disp('Gauss Hermite');
    S_GH = create_sparse_grid(N,w,knots_GH,@lev2knots_lin, @(i) sum(i-1), S_GH_old);
    Sr_GH = reduce_sparse_grid(S_GH);    
    evals_GH = evaluate_on_sparse_grid(f, S_GH, Sr_GH, evals_GH_old, S_GH_old, Sr_GH_old);
    err_GH(w) = max(abs(f_sampled - interpolate_on_sparse_grid(S_GH,Sr_GH,evals_GH,sampleset)));
    evals_GH_old = evals_GH;  
    S_GH_old = S_GH;
    Sr_GH_old = Sr_GH;       
    nb_pts_GH(w) = Sr_GH.size;
end


% repeat for GK knots

%knots_GK = @(n) knots_GK(n);   

w_max_GK=4;

knots_GK = @(n) knots_GK(n,0,1); 
err_GK=zeros(1,w_max_GK);
nb_pts_GK =zeros(1,w_max_GK);

S_GK_old=[];
Sr_GK_old=[];
evals_GK_old=[];



for w=1:w_max_GK       
    
    disp(w);
    disp('GK knots');
    S_GK = create_sparse_grid(N,w,knots_GK,@lev2knots_GK, @(i) sum(i-1), S_GK_old);
    Sr_GK = reduce_sparse_grid(S_GK);
    evals_GK = evaluate_on_sparse_grid(f, S_GK, Sr_GK, evals_GK_old, S_GK_old, Sr_GK_old);
    err_GK(w) = max(abs(f_sampled - interpolate_on_sparse_grid(S_GK,Sr_GK,evals_GK,sampleset)));   
    evals_GK_old = evals_GK;
    S_GK_old = S_GK;
    Sr_GK_old = Sr_GK;
    
    nb_pts_GK(w) = Sr_GK.size;
end

% convergence plots

figure
loglog(nb_pts_Lj, err_Lj,'-xr','LineWidth',2,'DisplayName','Normal Leja pts')
grid on
hold on
semilogy(nb_pts_GH, err_GH,'-ob','LineWidth',2,'DisplayName','Gauss Hermite pts')
semilogy(nb_pts_GK, err_GK,'-^k','LineWidth',2,'DisplayName','GK pts')

legend show
