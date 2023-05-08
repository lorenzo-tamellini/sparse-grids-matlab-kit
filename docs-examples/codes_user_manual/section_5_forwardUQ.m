
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The Sparse Grids Matlab kit user manual"
% C. Piazzola, L. Tamellini, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Forward UQ example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% spatial domain ----------------------------------- 
ab = [0,1];
N_el = 200; 
h = (ab(2)-ab(1))/N_el;
xx = linspace(ab(1),ab(2),N_el+1); 
% --------------------------------------------------

% for the definition of the random field -----------
N = 2; 
mu  = 1; 
s = [0.5, 0.1]; 
% --------------------------------------------------

PDE_rhs = @(x) ones(size(x)); 
I = @(y) QoI(xx,y,mu,s,PDE_rhs); 

% for the sparse grid ------------------------------

knots = @(n) knots_CC(n,-sqrt(3),sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 

%% loop over w 
w_vec = 2:8; 
exp_I = zeros(size(w_vec)); 
S_old = [];
Sr_old = []; 
evals_old = []; 
counter = 0; 
for w = [w_vec,10]
    counter = counter+1; 
    S = create_sparse_grid(N,w,knots,lev2knots,rule,S_old);
    Sr = reduce_sparse_grid(S);
    [exp_I(counter),evals] = quadrature_on_sparse_grid(I,S,Sr,evals_old,S_old,Sr_old);
    S_old = S; Sr_old = Sr;  evals_old = evals; 
end
err_exp_I = abs(exp_I(1:7)-exp_I(end)); 

% Fig. 11a 
figure
semilogy(w_vec,err_exp_I,'-*','LineWidth',2)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$w$','FontSize',18,'Interpreter','latex')
ylabel('Error','FontSize',18,'Interpreter','latex')

%% fix w = 4 

w = 4; 
S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

% computing the mean --------------------------------

[exp_I,I_on_Sr] = quadrature_on_sparse_grid(I,Sr); % takes in input the function

% computing the variance ----------------------------
int_I2 = quadrature_on_sparse_grid(I_on_Sr.^2,Sr); % int_{\Gamma} I^2 \rho
var_I = int_I2-exp_I^2;

% sparse grid approximation ------------------------- 

% Fig. 11b
figure
domain = [-sqrt(3),-sqrt(3);sqrt(3),sqrt(3)];
plot_sparse_grids_interpolant(S,Sr,domain,I_on_Sr)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')

% pdf ----------------------------------------------

% M samples in \Gamma
M = 5000; 
y_samples = -sqrt(3)+2*sqrt(3)*rand(N,M); 

% evaluations of the surrogate model
I_vals_rand = interpolate_on_sparse_grid(S,Sr,I_on_Sr,y_samples); 
[pdf_vals,pdf_pts] = ksdensity(I_vals_rand); % pdf

%  Fig. 11c
figure
H = histogram(I_vals_rand,'Normalization','pdf','LineWidth',1); 
hold on 
plot(pdf_pts,pdf_vals,'LineWidth',2)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

%% Sobol indices ----------------------------------
domain = [-sqrt(3) -sqrt(3); sqrt(3) sqrt(3)]; 
[Sob_idx, Tot_Sob_idx] = compute_sobol_indices_from_sparse_grid(S,Sr,I_on_Sr,domain,'legendre');
