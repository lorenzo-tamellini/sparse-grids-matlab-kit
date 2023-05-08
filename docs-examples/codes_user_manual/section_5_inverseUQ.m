
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The Sparse Grids Matlab kit user manual"
% C. Piazzola, L. Tamellini, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% Inverse UQ example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% spatial domain ----------------------------------- 
ab = [0,1];
% --------------------------------------------------

% for the definition of the random field -----------
N = 2; 
mu  = 1; 
s = [0.5, 0.5]; 
% --------------------------------------------------

PDE_rhs = @(x) ones(size(x)); 

%% noisy measurements for inversion 

% true solution ------------------------------------
K = 80;
x_k = linspace(ab(1),ab(2),K+2);
y_star = [0.9,-1.1];
u_star = PDE_solver(x_k,y_star,mu,s,PDE_rhs);
u_star = u_star(2:end-1); % select the values at the internal knots of the mesh
% --------------------------------------------------

% generate noisy data ------------------------------
sigma_eps = 0.01;
u_tilde = u_star + sigma_eps*randn(K,1);

% Fg. 12a
figure
plot(x_k(2:end-1),u_tilde,'*','LineWidth',2)
hold on 
plot(x_k,[0;u_star;0],'LineWidth',2)
legend('$\widetilde{u}_k$','$u(x,y^*)$','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

%% NLL 

% define the negative log-likelihood
u = @(y) PDE_solver(x_k,y,mu,s,PDE_rhs); 
misfits = @(y) [0;u_tilde;0] - u(y); % adding the boundary cond. to u_tilde just for simplicity
NLL = @(y) sum( misfits(y).^2  / (2*sigma_eps^2)) + (K-2)*log(sigma_eps^2) + (K-2)*0.5*log(2*pi);

y1_temp = linspace(0.2,sqrt(3),100); % interval containing y_star(1) = 0.9 
y2_temp = linspace(-0.5,-1.4,100); % interval containing y_star(2) = -1.1 
[Y1,Y2] = meshgrid(y1_temp,y2_temp);

yy_plot = [Y1(:)';Y2(:)']; % dim 2x10000 - each column is a evaluation point

NLL_plot = zeros(length(yy_plot),1);
for i = 1:length(yy_plot)
    NLL_plot(i) = NLL(yy_plot(:,i));
end

% Fig. 12b
figure
contourf(Y1,Y2,reshape(NLL_plot,100,100),100)
hold on
p = plot(y_star(1),y_star(2),'xr','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','r','DisplayName','$y^*$');
c = colorbar;
c.TickLabelInterpreter = 'latex';
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')
legend(p,'interpreter','latex')

% Fig. 12c
figure
surf(Y1,Y2,reshape(NLL_plot,100,100))
c = colorbar;
c.TickLabelInterpreter = 'latex';
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')

%% inversion 

% construct the surrogate model --------------------

knots = @(n) knots_CC(n,-sqrt(3), sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 

w = 5; 
S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

u_on_Sr = evaluate_on_sparse_grid(u,Sr); % u_on_Sr has dimension K x Sr.size:
                                         % the ith row contains u(x_i,Sr.knots)
u_approx = @(y) interpolate_on_sparse_grid(S,Sr,u_on_Sr,y); 
% note that interpolate_on_sparse_grids supports vector-valued functions,
% so in one call we get the values of u at all nodes x_i
LS = @(y) sum( ([0;u_tilde;0] - u_approx(y)).^2 );

% maximum a-posteriori estimate (MAP) -----------------------

% initial guess of the minimization in the center of \Gamma
y_start = [0;0]; 

% minimization of the NLL to find the MAP
y_MAP = fminsearch(@(y) LS(y),y_start); 

% estimate sigma_noise ------------------------------------
sigma_eps_approx = sqrt(mean(misfits(y_MAP).^2)); 

disp('========================')
disp('y_star')
disp(y_star)
disp('y_MAP')
disp(y_MAP')
disp('sigma_eps')
disp(sigma_eps)
disp('sigma_eps_approx') 
disp(sigma_eps_approx)

% covariance matrix --------------------------------------

% compute the Jacobian of the misfits wrt y
domain = [-sqrt(3), -sqrt(3); sqrt(3), sqrt(3)]; 
Jac_at_MAP = zeros(K-2,N); 
for i = 1:(K-2)
    Jac_at_MAP(i,:) = derive_sparse_grid(S,Sr,u_on_Sr(i,:),domain,y_MAP)'; 
end
Sigma_post = sigma_eps_approx^2*inv(Jac_at_MAP'*Jac_at_MAP);

% plot of posterior distribution -------------------------
x1_temp = linspace(y_MAP(1)-0.3,y_MAP(1)+0.3,100); 
x2_temp = linspace(y_MAP(2)-0.3,y_MAP(2)+0.3,100); 
[X1,X2] = meshgrid(x1_temp,x2_temp);
X = [X1(:) X2(:)];
pdf_post_approx = mvnpdf(X,y_MAP',Sigma_post);

% Fig. 13a
figure
contourf(X1,X2,reshape(pdf_post_approx,100,100),20);
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')

% Fig. 13b
figure
surf(X1,X2,reshape(pdf_post_approx,100,100));
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')
xlim([y_MAP(1)-0.3,y_MAP(1)+0.3])
ylim([y_MAP(2)-0.3,y_MAP(2)+0.3])

%% forward UQ based on the posterior

% spatial domain  
ab = [0,1];
N_el = 200; 
h = (ab(2)-ab(1))/N_el;
xx = linspace(ab(1),ab(2),N_el+1); 

% construct the surrogate model -----------------------
w = 4; 
knots = @(n) knots_normal(n,0,1);
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_lin; 

S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

%  change of variables: y = H'*z+y_MAP
H = chol(Sigma_post);
I = @(z) QoI(xx,H'*z+y_MAP,mu,s,PDE_rhs); 
I_on_Sr = evaluate_on_sparse_grid(I,Sr); 

% sampling the surrogate model -----------------------
M = 5000; 
z_samples = randn(N,M); 
I_vals_randn = interpolate_on_sparse_grid(S,Sr,I_on_Sr,z_samples); 
[pdf_vals,pdf_pts] = ksdensity(I_vals_randn); 

% Fig. 14b
figure
plot(pdf_pts,pdf_vals,'LineWidth',2)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

%% forward UQ based on the prior 

knots = @(n) knots_CC(n,-sqrt(3),sqrt(3));
rule = @(ii) sum(ii-1);
lev2knots = @lev2knots_doubling; 
S = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr = reduce_sparse_grid(S);

I = @(y) QoI(xx,y,mu,s,PDE_rhs); 
I_on_Sr = evaluate_on_sparse_grid(I,Sr); 

M = 5000; 
y_samples = -sqrt(3)+2*sqrt(3)*rand(N,M); 

% evaluations of the surrogate model
I_vals_rand = interpolate_on_sparse_grid(S,Sr,I_on_Sr,y_samples); 
[pdf_vals,pdf_pts] = ksdensity(I_vals_rand); % pdf

%  Fig. 14a
figure
plot(pdf_pts,pdf_vals,'LineWidth',2)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

