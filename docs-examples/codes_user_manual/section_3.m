
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The Sparse Grids Matlab kit user manual"
% C. Piazzola, L. Tamellini, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Listing 1: basic creation of a sparse grid
clear 

knots = @(n) knots_uniform(n,0,1);
lev2knots = @lev2knots_doubling; 
I = [1 1;
     1 2;
     2 1;
     3 1];
S = create_sparse_grid_multiidx_set(I,knots,lev2knots); 
Sr = reduce_sparse_grid(S); 

% plot of the sparse grid (Fig. 2d)
figure
plot_sparse_grid(S,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
axis square
xticks(0:0.2:1)
xlim([0 1])
ylim([0 1])

% plots of the tensor grids (Fig. 2a-c)
for i = 1:length(S)
    figure
    plot_sparse_grid(S(i),[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
    axis square
    xticks(0:0.2:1)
    set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
    xlim([0 1])
    ylim([0 1])
end

%% Listing 2: creation of a multi-index set 

N = 2;
rule = @(ii) sum(ii-1);
w = 4; 
base = 1;
I = multiidx_gen(N,rule,w,base); 

% plot of the multi-index set
figure
plot_multiidx_set(I,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

%% Example 2: multi-index sets

N = 2;
w = 4;
base = 1; 
g = [1,2];

% %% I_sum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multi-index sets
rule_iso = @(ii) sum(ii-1); % same as [~,rule_iso] = define_functions_for_rule('TD',N); 
I_iso = multiidx_gen(N,rule_iso,w,base);

rule_aniso = @(ii) sum(g(1:length(ii)).*(ii-1)); % same as [~,rule_aniso] = define_functions_for_rule('TD',g);
I_aniso = multiidx_gen(N,rule_aniso,w,base);

% corresponding sparse grids
knots = @(n) knots_leja(n,0,1,'sym_line');
lev2knots = @lev2knots_2step; 
S_iso = create_sparse_grid_multiidx_set(I_iso,knots,lev2knots); 
S_aniso = create_sparse_grid_multiidx_set(I_aniso,knots,lev2knots); 

% Fig. 5a
figure
plot_multiidx_set(I_iso,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

figure
plot_sparse_grid(S_iso,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
axis square
xlim([0 1])
ylim([0 1])

% Fig. 5b
figure
plot_multiidx_set(I_aniso,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

figure
plot_sparse_grid(S_aniso,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
axis square
xlim([0 1])
ylim([0 1])

% %% I_max %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% multi-index sets
rule_iso = @(ii) max(ii-1);
I_iso = multiidx_gen(N,rule_iso,w,base);
rule_aniso = @(ii) max(g(1:length(ii)).*(ii-1));
I_aniso = multiidx_gen(N,rule_aniso,w,base);

% corresponding sparse grids
knots = @(n) knots_leja(n,0,1,'sym_line');
lev2knots = @lev2knots_2step; 
S_iso = create_sparse_grid_multiidx_set(I_iso,knots,lev2knots); 
S_aniso = create_sparse_grid_multiidx_set(I_aniso,knots,lev2knots); 

% Fig. 5c
figure
plot_multiidx_set(I_iso,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

figure
plot_sparse_grid(S_iso,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 1])
ylim([0 1])

% Fig. 5d
figure
plot_multiidx_set(I_aniso,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

figure
plot_sparse_grid(S_aniso,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 1])
ylim([0 1])

% % I_prod %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multi-index sets
rule_iso = @(ii) prod(ii);
I_iso = multiidx_gen(N,rule_iso,w,base);
rule_aniso = @(ii) prod((ii).^g(1:length(ii)));
I_aniso = multiidx_gen(N,rule_aniso,w,base);

% corresponding sparse grids
knots = @(n) knots_leja(n,0,1,'sym_line');
lev2knots = @lev2knots_2step; 
S_iso = create_sparse_grid_multiidx_set(I_iso,knots,lev2knots); 
S_aniso = create_sparse_grid_multiidx_set(I_aniso,knots,lev2knots); 

% Fig. 5e
figure
plot_multiidx_set(I_iso,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

figure
plot_sparse_grid(S_iso,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 1])
ylim([0 1])

% Fig. 5f
figure
plot_multiidx_set(I_aniso,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:6)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 6])
ylim([0 6])

figure
plot_sparse_grid(S_aniso,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([0 1])
ylim([0 1])

%% Listing 3: basic creation of a sparse grid given a multi-index set defined by a rule

N = 2;
w = 3; 
knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling;
rule = @(ii) sum(ii-1);
% or by using the convenience function define_function_for_rule:
% [lev2nodes,idxset] = define_functions_for_rule('SM',N);
S = create_sparse_grid(N,w,knots,lev2knots,rule); 
Sr = reduce_sparse_grid(S); 

% Fig. 6b-h
for i = 1:length(S)
    figure
    plot_sparse_grid(S(i),[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
    axis square
    set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
    xlim([0 1])
    ylim([0 1])
end

% Fig. 6a
figure
plot_sparse_grid(S,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
xlim([0 1])
ylim([0 1])
xticks(0:0.2:1)
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

%% Listing 4: basic creation of an adaptive sparse grid 
clear 

f = @(y) exp(sum(y));
N = 2; 
knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling; 
controls = struct('nested',true); % each field of this struct specifies an optional argument
                                  % to control the algorithm
Ad = adapt_sparse_grid(f,N,knots,lev2knots,[],controls); 

%% 1st snippet in Sect. 3.4.3
clear 

I = [1 1;
     1 2;
     2 1;
     3 1];
knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling; 
S = create_sparse_grid_multiidx_set(I,knots,lev2knots); 
I_new = [I;
		 1 3;
		 2 2;
		 4 1]; 
I_new = sortrows(I_new); 
S_new = create_sparse_grid_multiidx_set(I_new,knots,lev2knots,S);

% plotting multi-iondex set and corresponding grid
figure
subplot(2,2,1)
plot_multiidx_set(I,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:4)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
xlim([0 5])
ylim([0 5])
title('multi-index set')

subplot(2,2,2)
plot_sparse_grid(S,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
xlim([0 1])
ylim([0 1])
xticks(0:0.2:1)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
title('sparse grid')

subplot(2,2,3)
plot_multiidx_set(I_new,'sk','LineWidth',2,'MarkerSize',22,'MarkerFaceColor','b')
axis square
xticks(0:5)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
xlim([0 5])
ylim([0 5])
title('new multi-index set')

subplot(2,2,4)
plot_sparse_grid(S_new,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',12);
axis square
xlim([0 1])
ylim([0 1])
xticks(0:0.2:1)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
title('new sparse grid')

%% 2nd snippet in Sect. 3.4.3
clear 

f = @(y) exp(sum(y));
N = 2; 
knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling; 
controls = struct('nested',true,'prof_tol',1e-5); 
Ad = adapt_sparse_grid(f,N,knots,lev2knots,[],controls);
controls_next = struct('nested',true); % default tolerance 1e-14
Next = adapt_sparse_grid(f,N,knots,lev2knots,Ad,controls_next);

% plotting multi-iondex set and corresponding grid
figure
subplot(2,2,1)
plot_multiidx_set(Ad.private.G,'sk','LineWidth',2,'MarkerSize',15,'MarkerFaceColor','b')
axis square
xticks(0:7)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
xlim([0 7])
ylim([0 7])
title('multi-index set')

subplot(2,2,2)
plot_sparse_grid(Ad.S,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',2); 
axis square
xlim([0 1])
ylim([0 1])
xticks(0:0.2:1)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
title('sparse grid')

subplot(2,2,3)
plot_multiidx_set(Next.private.G,'sk','LineWidth',2,'MarkerSize',15,'MarkerFaceColor','b')
axis square
xticks(0:7)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
xlim([0 7])
ylim([0 7])
title('new multi-index set')

subplot(2,2,4)
plot_sparse_grid(Next.S,[],'ok','LineWidth',3,'MarkerFaceColor','k','MarkerSize',2); 
axis square
xlim([0 1])
ylim([0 1])
xticks(0:0.2:1)
set(gca,'TickLabelInterpreter', 'latex','FontSize',14)
title('new sparse grid')

%% 3rd snippet in Sect. 3.4.3
clear 

knots = @(n) knots_CC(n,0,1);
lev2knots = @lev2knots_doubling; 
I_in = [1 1;
	    1 2;
	    2 1;
	    3 1];
S_in = create_sparse_grid_multiidx_set(I_in,knots,lev2knots); 
coeff_in = combination_technique(I_in); % cefficients of the combination tech. applied to I_in
new_idx = [4 1];
S_new = create_sparse_grid_add_multiidx(new_idx,S_in,I_in,coeff_in,knots,lev2knots);  

