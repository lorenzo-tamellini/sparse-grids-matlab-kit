
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The Sparse Grids Matlab kit - a Matlab implementation of sparse grids for 
% high-dimensional function approximation and uncertainty quantification"
% C. Piazzola, L. Tamellini, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example 2

clear 

N_vec = [2,4,6];
w_vec = 1:8; 
knots = @(n) knots_CC(n,0,1);

f = @(x) exp(sum(x)); 

% % without recycling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_eval = NaN*ones(length(N_vec),length(w_vec)); 
i = 0; 
for N = N_vec
    i = i+1; 
    
    [lev2knots,rule] = define_functions_for_rule('SM',N);
    j = 0; 
    for w = w_vec
        j = j+1; 
        S = create_sparse_grid(N,w,knots,lev2knots,rule); 
        Sr = reduce_sparse_grid(S); 
        f_evals = evaluate_on_sparse_grid(f,Sr); 
        nb_eval(i,j) = Sr.size;
    end
end

figure
cc = lines(length(N_vec));
pp = [];
for i = 1:length(N_vec)
    p = semilogy(w_vec,nb_eval(i,:),'-*','Color',cc(i,:),'LineWidth',2,'DisplayName',['$N$=',num2str(N_vec(i))]);
    pp = [pp,p];
    hold on
end

% % with recycling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nb_eval_red = NaN*ones(length(N_vec),length(w_vec));
i = 0; 
for N = N_vec
    i = i+1; 
    
    [lev2knots,rule] = define_functions_for_rule('SM',N);
    S_old = []; % initialize the grids 
    Sr_old = [];
    evals_old = []; 

    j = 0;
    tic
    for w = w_vec
        j = j+1; 
        
        S = create_sparse_grid(N,w,knots,lev2knots,rule,S_old);  % recycle from S_old
        Sr = reduce_sparse_grid(S); 
        
        f_evals = evaluate_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old); 
        % [int_val, evals] = quadrature_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old); 
        % f_int(i,j) = int_val; 
        nb_eval_red(i,j) = length(f_evals)-length(evals_old);
        
        S_old = S; 
        Sr_old = Sr; 
        evals_old = f_evals; 
    end
    
end

for i = 1:length(N_vec)
    semilogy(w_vec,nb_eval_red(i,:),'--*','Color',cc(i,:),'LineWidth',2)
    hold on
end

grid on
legend(pp,'FontSize',16,'Location','northwest','Interpreter','latex')
xlabel('w','FontSize',18,'Interpreter','latex')
ylabel('Number of evaluations','Fontsize',18,'Interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlim([w_vec(1) w_vec(end)])
xticks(w_vec)

% plotting the percentage of saved computations
perc_rid = (nb_eval-nb_eval_red)./nb_eval*100;
figure
for i = 1:length(N_vec)
    plot(w_vec,perc_rid(i,:),'-*','Color',cc(i,:),'LineWidth',2,'DisplayName',['$N$=',num2str(N_vec(i))])
    hold on
end
grid on
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('w','FontSize',18,'Interpreter','latex')
ylabel('$\%$ of saved evaluations','FontSize',18,'Interpreter','latex')
legend('location','southeast','interpreter','latex')
xlim([w_vec(1),w_vec(end)])
xticks(w_vec)


