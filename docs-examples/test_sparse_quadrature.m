
%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


%% test 1, exact solution available

clear

% function to be integrate, in (-1,1)^N. input column points, out row vector

f = @(x,b) prod(1./sqrt(x+b));
b=3;
N = 6;
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b));
I_ex = I_1d^N;


% define sparse grid
[lev2knots,idxset]=define_functions_for_rule('TD',N);
knots=@(n) knots_uniform(n,-1,1,'nonprob');


% for loop

w_max=6;
q_error=zeros(1,w_max);
work=zeros(1,w_max);

S_old=[];
Sr_old=[];
evals_old=[];

tic
for w=0:w_max

    disp(w)

    % create grid
    [S,C]=create_sparse_grid(N,w,knots,lev2knots,idxset);
    Sr=reduce_sparse_grid(S);

    [I,evals_old]=quadrature_on_sparse_grid(@(x)f(x,b),S,Sr,evals_old,S_old,Sr_old);
    
    S_old=S;
    Sr_old=Sr;    
    
    %compute convergence error
    q_error(w+1)=abs(I_ex-I);

    % compute work
    work(w+1)=Sr.size;
        
end
toc

% error w.r.t. level w
figure
semilogy(0:w_max,q_error,'-or','DisplayName','Numerical Error, w.r.t. grid level');
hold on
semilogy(0:w_max,1./exp(0:w_max),'--o','DisplayName','exp(-level)')
legend show

% error w.r.t. nb. points
figure
loglog(work,q_error,'-or','DisplayName','Numerical Error, w.r.t. #points');
legend show

%% test 2, no exact solution available. Also, we take this chance and show an example with exponential random variables


clear

% dimension of space
N=2;

% we use a simple TD rule, up to this level 
w_max=15;

% function to be integrated

% the function resembles a discontinuous (step) function the higher the factor in the argument of the exp is
f = @(x) 1/(1+exp(0.7*sum(x))); 

% oscillatory function: change the factor to see different speed of convergence 
% the higher the factor the more oscillatory the function (and slower the
% convergence)
% f = @(x) cos(0.2*sum(x)); 

% peak function: the factor in front of the norm acts on the steepness of the peak 
% the higher the number the steeper the peak (and slower the convergence)
% f = @(x) 1/(1+0.1*norm(x)^2); 

lambda = 1; 
points = @(n) knots_exponential(n,lambda);

quad = zeros(1,w_max);
nb_pts = zeros(1,w_max);

% we introduce auxiliary containers to recycle previous evaluations and speed up the computation
S_old=[];
Sr_old=[];
evals_old=[];

% convergence loop
for w=1:w_max       
    
    S = create_sparse_grid(N,w,points,@lev2knots_lin, @(i) sum(i-1), S_old);
    Sr = reduce_sparse_grid(S);    
    [res, evals]  = quadrature_on_sparse_grid(f, S, Sr, evals_old, S_old, Sr_old);
    quad(w) = res;
    evals_old = evals;
    S_old = S;
    Sr_old = Sr;    
    nb_pts(w) = Sr.size;
    
end

% exact integral
S = create_sparse_grid(N,w_max+4,points,@lev2knots_lin,@(i) sum(i-1),S_old);
Sr = reduce_sparse_grid(S);
exact = quadrature_on_sparse_grid(f,S,Sr,evals_old,S_old,Sr_old);

% errors
err = abs(quad - exact);

figure
loglog(nb_pts, err,'-ob','LineWidth',2)
xlabel('nb. of quadrature points')
ylabel('error')
grid on



