
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The Sparse Grids Matlab kit - a Matlab implementation of sparse grids for 
% high-dimensional function approximation and uncertainty quantification"
% C. Piazzola, L. Tamellini, 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Computational cost - vary N 

clear 

N_vec = 2:10; 
w = 3; 
knots = @(n) knots_CC(n,0,1);

time = zeros(10,length(N_vec)); 
time_tot = zeros(10,length(N_vec)); 
sg_dim = zeros(size(N_vec));

for k = 1:10 % repeat 10 times
    c = 0; 
    for N = N_vec
        c = c+1; 
        tic 
        [lev2nodes,idxset] = define_functions_for_rule('SM',N);
        S = create_sparse_grid(N,w,knots,lev2nodes,idxset); 
        time(k,c) = toc; 
        Sr = reduce_sparse_grid(S);
        time_tot(k,c) = toc; 
        sg_dim(c) = Sr.size;
    end
end

time_avg = mean(time,1);
time_tot_avg = mean(time_tot,1);

w = 5; 
time_1 = zeros(10,length(N_vec)); 
time_1_tot = zeros(10,length(N_vec)); 
sg_dim_1 = zeros(size(N_vec));


for k = 1:10 % repeat 10 times
    c = 0; 
    for N = N_vec
        c = c+1; 
        tic
        [lev2nodes,idxset] = define_functions_for_rule('SM',N);
        S = create_sparse_grid(N,w,knots,lev2nodes,idxset);
        time_1(k,c) = toc; 
        Sr = reduce_sparse_grid(S);
        time_1_tot(k,c) = toc; 
        sg_dim_1(c) = Sr.size;
    end
end

time_1_avg = mean(time_1,1);
time_1_tot_avg = mean(time_1_tot,1);


cc = lines(2); 
figure
yyaxis left
semilogy(N_vec,time_tot_avg,'-*','Color',cc(1,:),'LineWidth',2)
hold on 
semilogy(N_vec,time_1_tot_avg,'-.s','Color',cc(1,:),'LineWidth',2)
grid on 
yticks(10.^(-3:2:1))
xlabel('N','FontSize',18,'Interpreter','latex')
ylabel('Computational time [s]','FontSize',18,'Interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
set(gca,'YMinorGrid','off')
hold off

yyaxis right
semilogy(N_vec,sg_dim,'-*','Color',cc(2,:),'LineWidth',2)
hold on
semilogy(N_vec,sg_dim_1,'-.s','Color',cc(2,:),'LineWidth',2)
p1 = semilogy(nan,nan,'-*k','LineWidth',2);
p2 = semilogy(nan,nan,'-.sk','LineWidth',2);
yticks(10.^(1:2:5))
ylabel('Size of the sparse grid','FontSize',18,'Interpreter','latex')
xlabel('N','FontSize',18,'Interpreter','latex')
legend([p1,p2],'$w=3$','$w=5$','interpreter','latex','fontsize',16,'Location','northwest'); 
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)

figure
plot(N_vec,(time_tot_avg-time_avg)./time_tot_avg*100,'-*','Color',cc(1,:),'LineWidth',2)
hold on 
plot(N_vec,(time_1_tot_avg-time_1_avg)./time_1_tot_avg*100,'-.s','Color',cc(1,:),'LineWidth',2)
grid on
legend('$w=3$','$w=5$','interpreter','latex','fontsize',16,'Location','northwest'); 
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('N','FontSize',18,'Interpreter','latex')
ylabel('$\%$ of time for reduction','FontSize',18,'Interpreter','latex')

%% Computational cost - vary w
  
clear 

w_vec = 2:10; 
N = 3; 
knots = @(n) knots_CC(n,0,1);
time = zeros(10,length(w_vec)); 
time_tot = zeros(10,length(w_vec)); 
sg_dim = zeros(size(w_vec));

for k = 1:10
    c = 0; 
    for w = w_vec
        c = c+1; 
        tic
        [lev2nodes,idxset] = define_functions_for_rule('SM',N);
        S = create_sparse_grid(N,w,knots,lev2nodes,idxset); 
        time(k,c) = toc; 
        Sr = reduce_sparse_grid(S);
        time_tot(k,c) = toc; 
        sg_dim(c) = Sr.size;
    end
end

time_avg = mean(time,1);
time_tot_avg = mean(time_tot,1);

N = 5; 
time_1 = zeros(10,length(w_vec)); 
time_1_tot = zeros(10,length(w_vec)); 
sg_dim_1 = zeros(size(w_vec));

for k = 1:10
    c = 0; 
    for w = w_vec
        c = c+1; 
        tic
        [lev2nodes,idxset] = define_functions_for_rule('SM',N);
        S = create_sparse_grid(N,w,knots,lev2nodes,idxset); 
        time_1(k,c) = toc; 
        Sr = reduce_sparse_grid(S);
        time_1_tot(k,c) = toc; 
        sg_dim_1(c) = Sr.size;
    end
end

time_1_avg = mean(time_1,1);
time_1_tot_avg = mean(time_1_tot,1);

cc = lines(2);
figure
yyaxis left
semilogy(w_vec,time_tot_avg,'-*','Color',cc(1,:),'LineWidth',2)
hold on
semilogy(w_vec,time_1_tot_avg,'-.s','Color',cc(1,:),'LineWidth',2)
grid on 
yticks(10.^(-3:2:1))
ylim([10^(-3) 10^1])
ylabel('Computational time [s]','FontSize',18,'Interpreter','latex')
xlabel('w','FontSize',18,'Interpreter','latex')
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
set(gca,'YMinorGrid','off')
hold off

yyaxis right
semilogy(w_vec,sg_dim,'-*','Color',cc(2,:),'LineWidth',2)
hold on 
semilogy(w_vec,sg_dim_1,'-.s','Color',cc(2,:),'LineWidth',2)
grid on 
ylabel('Size of the sparse grid','FontSize',18,'Interpreter','latex')
xlabel('w','FontSize',18,'Interpreter','latex')
p1 = plot(nan,nan,'-*k','LineWidth',2);
p2 = plot(nan,nan,'-.sk','LineWidth',2);
legend([p1,p2],'$N=3$','$N=5$','interpreter','latex','fontsize',16,'location','northwest'); 
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
set(gca,'YMinorGrid','off')
yticks(10.^(0:2:6))

figure
plot(w_vec,(time_tot_avg-time_avg)./time_tot_avg*100,'-*','Color',cc(1,:),'LineWidth',2)
hold on 
plot(w_vec,(time_1_tot_avg-time_1_avg)./time_1_tot_avg*100,'-.s','Color',cc(1,:),'LineWidth',2)
grid on
legend('$N=3$','$N=5$','interpreter','latex','fontsize',16,'Location','northwest'); 
set(gca,'TickLabelInterpreter', 'latex','FontSize',18)
xlabel('w','FontSize',18,'Interpreter','latex')
ylabel('$\%$ of time for reduction','FontSize',18,'Interpreter','latex')

