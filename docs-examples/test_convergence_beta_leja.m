%-------------------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2018 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%-------------------------------------------------------------

close all

%% Plot the weights of each Beta-Leja quadrature rule 
% Beta distribution (beta=1): rho(x)=x^alpha*(1-x)^beta 

clear

alpha = 0.3; 
beta = -0.2; 

[~,W]=compute_beta_leja_knots_and_weights_50(alpha,beta);

figure
semilogy(abs(W),'-') % plots weights of each quadrature rule

figure
plot(W','-') % plots trend of weights of each quadrature rule


%% Testing each Beta-Leja quadrature on computation of moments of a Beta random variable

clear

alpha = 5; 
beta = -0.5;

% Moments of a Beta distribution rho(x)=x^alpha*(1-x)^beta: 
% gamma(alpha+1+p)*gamma(alpha+beta+2)/(gamma(alpha+1)*gamma(alpha+beta+p+2))
p_max = 50;
mom = zeros(1,1+p_max);
for p = 0:1:p_max
    mom(p+1) = gamma(alpha+1+p)*gamma(alpha+beta+2)/(gamma(alpha+1)*gamma(alpha+beta+p+2));
end

% Quadrature error for polynomials
imax = 50;
err = zeros(1+p_max,imax);
errGJac = zeros(1+p_max,imax);

% Beta-Leja knots and weights
[X,W] = compute_beta_leja_knots_and_weights_50(alpha,beta); 


% for each formula, we test its approximation of increasing moments
for n=1:50
    % Beta-Leja quadrature rule using n nodes
    % select the first n knots and the corresponding weights from X,W
    % compute above
    [x_Lj,w_Lj]=knots_general_weighted_leja(n,X,W);
        
    % Gauss-Jacobi quadrature of same accuracy
    [x_GJac,w_GJac] = knots_beta(ceil(n/2),alpha,beta); 
    
    for p=0:p_max
        if p<n+5 % if the degree is "not too much" compute error
            err(1+p,n) = abs(mom(1+p) - dot(x_Lj.^p,w_Lj) );
            errGJac(1+p,n) = abs(mom(1+p) - dot(x_GJac.^p,w_GJac) );
        else % otherwise,  error is just too much,  we  set it to NaN
            err(1+p,n) = NaN;
            errGJac(1+p,n) = NaN;    
        end
    end
    
    % Plotting errors
    figure
    semilogy(0:p_max, err(:,n),'o',0:p_max, errGJac(:,n),'+','LineWidth',2)
    grid on
    set(gca,'FontSize',14)
    legend(sprintf('Beta-Leja (n=%i)',n),sprintf('Gauss-Jacobi (n=%i)',ceil(n/2)))
    set(legend,'Location','NorthWest')
    xlabel('Polynomial degree p')
    ylabel('Absolute quadrature error')
    
    pause
    
end

%% 1d quadrature - convergence test: increase number of points in the quadrature rule

clear

alpha = 1; 
beta = 3; 

imax=50;

% function to be integrated

% the function resembles a discontinuous (step function) the higher the factor in the argument of the exp is
f = @(x) 1./(1+exp(0.5*x)); 

% oscillatory function: change the period to see different speed of convergence 
% the higher the factor the more oscillatory the function 
% f = @(x) cos(0.5*x); 

% peak function: the factor in front of the norm acts on the steepness of the peak 
% the higher the number the steeper the peak
% f = @(x) 1./(1+0.1*x.^2); 

quad_Lj = zeros(1,imax);
quad_GJac = zeros(1,imax);
nb_pts =zeros(1,imax);

% Leja knots 
[X,W]=compute_beta_leja_knots_and_weights_50(alpha,beta);

% refining quad rule
for i=1:imax
    
    n = lev2knots_lin(i);
    nb_pts(i) = n;

    [x_Lj,w_Lj] = knots_general_weighted_leja(n,X,W);   
    [x_GJac,w_GJac] = knots_beta(n,alpha,beta);
    
    quad_Lj(i) = dot(f(x_Lj),w_Lj);
    quad_GJac(i) = dot(f(x_GJac),w_GJac);
end

% exact integral
[x_GJac,w_GJac] = knots_beta(100,alpha,beta);
exact = dot(f(x_GJac),w_GJac);
err_Lj = abs(quad_Lj - exact);
err_GJac = abs(quad_GJac - exact);


% Plotting errors
figure
semilogy(nb_pts, err_Lj,'-xr','LineWidth',2,'DisplayName','Exponential-Leja pts')
grid on
hold on
semilogy(nb_pts, err_GJac,'-ob','LineWidth',2,'DisplayName','Gauss-Jacobi pts')

legend show
set(legend,'Location','SouthWest')


ylim([1e-16 10])


%% 2d quadrature - convergence test

clear

alpha = 3;
beta = 1;

% dimension of space
N=2;

% we use a simple TD rule, up to this level
w_max=15;

% function to be integrated
f = @(x) 1/(1+exp(0.1*sum(x))); 
% f = @(x) cos(0.2*sum(x)); 
% f = @(x) 1/(1+0.1*norm(x)^2); 

% Leja knots 
[X,W]=compute_beta_leja_knots_and_weights_50(alpha,beta);

knots_Lj = @(n) knots_general_weighted_leja(n,X,W);   
knots_GJac = @(n) knots_beta(n,alpha,beta);

quad_Lj=zeros(1,w_max);
quad_GJac=zeros(1,w_max);

nb_pts_Lj =zeros(1,w_max);
nb_pts_GJac =zeros(1,w_max);

% we introduce auxiliary containers to recycle previous evaluations and speed up the computation
S_Lj_old=[];
Sr_Lj_old=[];
evals_Lj_old=[];

S_GJac_old=[];
Sr_GJac_old=[];
evals_GJac_old=[];

% the convergence loop for Leja and Gauss-Jacobi
for w=1:w_max       
    
    disp('Beta-Leja');
    S_Lj = smolyak_grid(N,w,knots_Lj,@lev2knots_2step, @(i) sum(i-1), S_Lj_old); % using 2step rule to ramp faster to rules with high number of points
    Sr_Lj = reduce_sparse_grid(S_Lj);
    [res,evals]= quadrature_on_sparse_grid(f, S_Lj, Sr_Lj, evals_Lj_old, S_Lj_old, Sr_Lj_old);
    quad_Lj(w) = res;
    evals_Lj_old = evals;
    S_Lj_old=S_Lj;
    Sr_Lj_old = Sr_Lj;
    nb_pts_Lj(w) = Sr_Lj.size;

    disp('Gauss-Jacobi');
    S_GJac = smolyak_grid(N,w,knots_GJac,@lev2knots_lin, @(i) sum(i-1), S_GJac_old);
    Sr_GJac = reduce_sparse_grid(S_GJac);    
    [res, evals]  = quadrature_on_sparse_grid(f, S_GJac, Sr_GJac, evals_GJac_old, S_GJac_old, Sr_GJac_old);
    quad_GJac(w) = res;
    evals_GJac_old = evals;
    S_GJac_old = S_GJac;
    Sr_GJac_old = Sr_GJac;    
    nb_pts_GJac(w) = Sr_GJac.size;
end


% exact integral
disp('computing reference solution');
S_GJac = smolyak_grid(N,w_max+4,knots_GJac,@lev2knots_lin, @(i) sum(i-1), S_GJac_old);
Sr_GJac = reduce_sparse_grid(S_GJac);
exact = quadrature_on_sparse_grid(f,S_GJac,Sr_GJac,evals_GJac_old,S_GJac_old,Sr_GJac_old);

% errors
err_Lj=abs(quad_Lj - exact);
err_GJac=abs(quad_GJac - exact);

figure
loglog(nb_pts_Lj, err_Lj,'-xr','LineWidth',2,'DisplayName','Beta-Leja pts')
grid on
hold on
loglog(nb_pts_GJac, err_GJac,'-ob','LineWidth',2,'DisplayName','Gauss Jacobi pts')

legend show
set(legend,'Location','SouthWest')

%% interpolation error in 1D

clear

alpha = -0.2; 
beta = -0.8; 

% Leja knots 
[X,W]=compute_beta_leja_knots_and_weights_50(alpha,beta);

% function to be interpoled
f1 = @(x) 1./(1+0.1*x.^2); 
f2 = @(x) 1./(2+exp(x)); 
f3 = @(x) cos(3*x+1); 

ff = {f1;f2;f3}; 

figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])

for j=1:3
    
    f = ff{j}; 

    % sample the function 
    samplesize = 50;
    sampleset = (linspace(0,1,samplesize)); 
    F_samples = f(sampleset);

    % the convergence analysis
    imax=50;

    err_Lj=zeros(1,imax);
    err_GH=zeros(1,imax);

    nb_pts =zeros(1,imax);

    for i=1:imax

        n = lev2knots_lin(i);
        nb_pts(i) = n;
        nnn = 1:n;

        % here we build the lagrange interpolant for Leja and evaluate the error
        [x_Lj,w_Lj] = knots_general_weighted_leja(n,X,W);       
        interp_Lj = zeros(1,samplesize);
        for k=nnn
            interp_Lj =  interp_Lj + f(x_Lj(k))*lagr_eval(x_Lj(k), x_Lj(nnn~=k),sampleset);
        end    
        err_Lj(i) = max(abs(F_samples - interp_Lj)); 

        % repeat for Gauss Hermite
        [x_GJac,w_GJac]=knots_beta(n,alpha,beta);
        interp_GJac= zeros(1,samplesize);
        for k=nnn
            interp_GJac =  interp_GJac + f(x_GJac(k))*lagr_eval(x_GJac(k), x_GJac(nnn~=k),sampleset);
        end    
        err_GJac(i) = max(abs(F_samples - interp_GJac)); 


    end
    
    subplot(2,3,j)
    plot(sampleset, F_samples,'-k','LineWidth',1.5,'DisplayName','f')
    hold on 
    plot(sampleset, interp_Lj,'*r','DisplayName','interp - Leja pts ')
    plot(sampleset, interp_GJac,'+b','DisplayName','interp - Gauss pts ')
    legend show
    set(legend,'Location','SouthWest')
    title(strrep(char(ff{j}),'@(x)','f='))

    % Plotting errors
    subplot(2,3,3+j)
    semilogy(nb_pts, err_Lj,'-xr','LineWidth',1.5,'DisplayName','Beta-Leja pts')
    hold on
    semilogy(nb_pts, err_GJac,'-ob','LineWidth',1.5,'DisplayName','Gauss-Jacobi pts')
    grid on
    legend show
    set(legend,'Location','SouthWest')

end
