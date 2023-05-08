
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The Sparse Grids Matlab kit user manual"
% C. Piazzola, L. Tamellini, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% First 4 snippets of Sect. 4.1
clear 

N = 2;
w = 3; 
knots = @(n) knots_CC(n,0,1);
[lev2knots,rule] = define_functions_for_rule('SM',N);
S = create_sparse_grid(N,w,knots,lev2knots,rule); 
Sr = reduce_sparse_grid(S);

% %% 1st snippet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = @(x) exp(sum(x)); 
f_on_Sr = evaluate_on_sparse_grid(f,Sr); 

% Fig. 7a
figure
plot3(Sr.knots(1,:),Sr.knots(2,:),f_on_Sr,'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',8)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')

% %% 2nd snippet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = quadrature_on_sparse_grid(f_on_Sr,Sr); 
q2 = quadrature_on_sparse_grid(f_on_Sr.^2,Sr);

% %% 3rd snippet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q,f_on_Sr] = quadrature_on_sparse_grid(f,Sr); 

% %% 4th snippet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% building the Cartesian grid of 15x15 equispaced points
y = linspace(0,1,15); 
[Y1,Y2] = meshgrid(y,y);
% rearrange the points such that each point is a column
eval_points = [Y1(:)';Y2(:)']; % matrix of dim. 2x15
% S, Sr and f_on_Sr as above  
f_vals = interpolate_on_sparse_grid(S,Sr,f_on_Sr,eval_points);

% Fig. 7b
figure
surf(Y1,Y2,reshape(f_vals,15,15))
hold on 
plot3(eval_points(1,:),eval_points(2,:),f_vals,'o','Color','r','LineWidth',1,'MarkerFaceColor','r','MarkerSize',5)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')

%% Last snippet of Sect. 4.1

clear 

w = 3; 
knots = @(n) knots_CC(n,0,1);

% N = 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
N = 3; 
f = @(x) exp(sum(x));

[lev2knots,rule] = define_functions_for_rule('SM',N);
S = create_sparse_grid(N,w,knots,lev2knots,rule); 
Sr = reduce_sparse_grid(S); 

% Fig. 8a
figure
plot3_sparse_grid(S,[],'ok','MarkerFaceColor','k','MarkerSize',8);
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')
zlabel('$y_3$','FontSize',18,'Interpreter','latex')

f_on_Sr = evaluate_on_sparse_grid(f,Sr); 

domain = [0 0 0; 1 1 1];

% Fig. 8b
figure
plot_sparse_grids_interpolant(S,Sr,domain,f_on_Sr,'nb_contourfs',5,'nb_countourf_lines',10)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('$y_1$','FontSize',18,'Interpreter','latex')
ylabel('$y_2$','FontSize',18,'Interpreter','latex')
zlabel('$y_3$','FontSize',18,'Interpreter','latex')
c = colorbar;
c.TickLabelInterpreter = 'latex';

% N = 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 4;
f = @(x) exp(x(1,:)+0.5*x(2,:)+1.5*x(3,:)+2*x(4,:));
domain = [0 0 0 0; 1 1 1 1]; 

[lev2knots,rule] = define_functions_for_rule('SM',N);
S = create_sparse_grid(N,w,knots,lev2knots,rule); 
Sr = reduce_sparse_grid(S); 

f_on_Sr = evaluate_on_sparse_grid(f,Sr); 

% Fig. 9
figure
plot_sparse_grids_interpolant(S,Sr,domain,f_on_Sr,'two_dim_cuts',[1 2 3 4 1 4])

%% Snippet of Sect. 4.2.1

clear 

N = 2;
f = @(x) exp(sum(x)); 
knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling; 
rule = @(ii) sum(ii-1);
% consider a sequence of grids with increasing level
w_vec = 1:6; 
S_old = []; % initialize the grid
Sr_old = []; % initialize the reduced grid
evals_old = []; 
for w = w_vec
    % recycle tensor grids from _old        
    S = create_sparse_grid(N,w,knots,lev2knots,rule,S_old);  
    Sr = reduce_sparse_grid(S); 
    % recycle evaluations on S_old        
    f_evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old); 
    % recycling with quadrature would be obtained in the same way
    % [int_val, evals] = quadrature_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old); 
    % 
    % do something with f_evals
    % 
    % then update the containers
    S_old = S; 
    Sr_old = Sr; 
    evals_old = f_evals; 
end


%% Snippet Sect. 4.4 - Convert from sparse grid to PCE
clear 

f = @(y) exp(sum(y)); 

N = 2;

% % generate the sparse grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 3;
base = 1; 
knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling;
rule = @(ii) sum(ii-1);
I = multiidx_gen(N,rule,w,base);
S = create_sparse_grid_multiidx_set(I,knots,lev2knots);
Sr = reduce_sparse_grid(S); 

% % evaluations of the function f at the grid points %%%%%%%%%%%%%%%%%%%%%%
f_on_Sr = evaluate_on_sparse_grid(f,Sr); 

% % compute gPCE coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
domain = [0 0; 1 1]; 
[PCE_coeffs,Lambda] = convert_to_modal(S,Sr,f_on_Sr,domain,'legendre'); 

% evaluate the gPCE
eval_points = rand(N,10); % points where to evaluate the gPCE of f
PCE_vals = zeros(1,length(eval_points)); 
for i = 1:length(PCE_coeffs)
    PCE_vals = PCE_vals+PCE_coeffs(i)*lege_eval_multidim(eval_points,Lambda(i,:),0,1); 
end

% Fig. 10a 
figure
plot_multiidx_set(I,'sk','LineWidth',2,'MarkerSize',15,'MarkerFaceColor','b')
axis square
xticks(0:9)
yticks(0:9)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 9])
ylim([0 9])

% fig. 10b
figure
plot_multiidx_set(Lambda,'sk','LineWidth',2,'MarkerSize',15,'MarkerFaceColor','b')
axis square
xticks(0:9)
yticks(0:9)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 9])
ylim([0 9])
