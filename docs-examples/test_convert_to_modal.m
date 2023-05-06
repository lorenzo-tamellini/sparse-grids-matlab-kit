%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2022 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


%% test convert to modal, using HC grid in 2D, interpolating Legendre polynomial

clear

% the sparse grid
N=2; 
w=5; 
knots=@(n) knots_uniform(n,-1,1,'nonprob'); 
lev2knots=@lev2knots_lin; 
idxset=@(i) prod(i); 

S=create_sparse_grid(N,w,knots,lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% the domain of the grid
domain=[-ones(1,N); ones(1,N)];


% compute a legendre polynomial over the sparse grid
X=Sr.knots;
nodal_values = 4*lege_eval_multidim(X,[4 0],-1,1)+2*lege_eval_multidim(X,[1 1],-1,1);

% conversion from the points to the legendre polynomial. I should recover it exactly
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'legendre');

[K,modal_coeffs]

%% test convert_to_modal, here each random var is defined over a different interval

clear

% left-ends of the intervals
a=[0 -1];
% right-ends of the intervals
b=[2 1];


% sparse grid
N=2; w=5; 

knots1=@(n) knots_uniform(n,a(1),b(1),'nonprob'); 
knots2=@(n) knots_uniform(n,a(2),b(2),'nonprob'); 
lev2knots=@lev2knots_lin; 
idxset=@(i) sum(i-1);

[S,C]=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% evaluate a linear combination of Legendre polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*lege_eval_multidim(X,[3 0],a,b)+7*lege_eval_multidim(X,[1 1],a,b);


domain=[a;b];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'legendre');

[K,modal_coeffs]


%% Hermite case
clear

flag='hermite';

% means
mu=[1 -1];
% standard deviations
sig=[2 3];

% The sparse grid. Note that here I don't really care what points I have used to build the nodal
% interpolant. They could even be evenly spaced points! It's just a bunch of evaluations of a function.
% of course bad points make the matrix ill conditioned
N=2; w=5; 

knots1=@(n) knots_normal(n,mu(1),sig(1));
knots2=@(n) knots_normal(n,mu(2),sig(2));
lev2knots=@lev2knots_lin; 
idxset=@(i) sum(i-1);

S=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% Same procedure as before, now with a linear comb of Hermite polynomials

X=Sr.knots;

nodal_values = 6-4*herm_eval_multidim(X,[3 1],mu,sig)+7*herm_eval_multidim(X,[1 2],mu,sig);


domain=[mu;sig];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,flag);

[K,modal_coeffs]


%% chebyshev case

clear

% left-ends of the intervals
a=[0 -1];
% right-ends of the intervals
b=[2 1];


% sparse grid
N=2; 
w=5; 
knots1=@(n) knots_uniform(n,a(1),b(1),'nonprob'); 
knots2=@(n) knots_uniform(n,a(2),b(2),'nonprob');
lev2knots=@lev2knots_lin; idxset=@(i) sum(i-1);

[S,C]=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% evaluate a linear comb of Chebyshev polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*cheb_eval_multidim(X,[3 2],a,b)+ 7*cheb_eval_multidim(X,[1 1],a,b);


domain=[a;b];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'chebyshev');

[K,modal_coeffs]

%% laguerre case

clear

% parameter
lambda = [1,2]; 

% sparse grid
N=2; 
w=5; 
knots1=@(n) knots_exponential(n,lambda(1)); 
knots2=@(n) knots_exponential(n,lambda(2));
lev2knots=@lev2knots_lin; idxset=@(i) sum(i-1);

[S,C]=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% evaluate a linear comb of Chebyshev polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*lagu_eval_multidim(X,[3 2],lambda)+ 7*lagu_eval_multidim(X,[2 1],lambda);


domain=lambda;
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'laguerre');

[K,modal_coeffs]


%% generalized laguerre case

clear

% parameter
alpha = [1 2]; 
beta = [3 4]; 

% sparse grid
N=2; 
w=5; 
knots1=@(n) knots_gamma(n,alpha(1),beta(1)); 
knots2=@(n) knots_gamma(n,alpha(2),beta(2));
lev2knots=@lev2knots_lin; idxset=@(i) sum(i-1);

[S,C]=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% evaluate a linear comb of Chebyshev polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*generalized_lagu_eval_multidim(X,[3 2],alpha,beta)+ 7*generalized_lagu_eval_multidim(X,[4 1],alpha,beta);

domain=[alpha;beta];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'generalized laguerre');

[K,modal_coeffs]

%% Jacobi case
clear

% parameter
alpha = [1 1]; 
beta = [3 4]; 
a = [-1 -1]; 
b = [1 1];

% sparse grid
N=2; 
w=5; 
knots1=@(n) knots_beta(n,alpha(1),beta(1),a(1),b(1)); 
knots2=@(n) knots_beta(n,alpha(2),beta(2),a(2),b(2));
lev2knots=@lev2knots_lin; idxset=@(i) sum(i-1);

[S,C]=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% evaluate a linear comb of Chebyshev polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*jacobi_prob_eval_multidim(X,[3 2],alpha,beta,a,b)+ 7*jacobi_prob_eval_multidim(X,[4 1],alpha,beta,a,b);

domain=[alpha;beta;a;b];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'jacobi');

[K,modal_coeffs]

%% test using different polynomials in different directions

clear
clc

a=0;
b=1.5;
mu = 1;
sig = 0.2;
a2=1.4;
b2=3.2;
lambda=1.5;  

N=4; w=7; 
knots={ @(n) knots_uniform(n,a,b,'prob'), @(n) knots_normal(n,mu,sig),...
        @(n) knots_uniform(n,a2,b2,'prob'), @(n) knots_exponential(n,lambda)};
lev2knots={@lev2knots_lin, @lev2knots_lin, @lev2knots_lin, @lev2knots_lin}; 
idxset=@(i) sum(i-1);
[S,C]=create_sparse_grid(N,w,knots,lev2knots,idxset);
Sr=reduce_sparse_grid(S);



X=Sr.knots;

ll1=[2 3 1 1];
ll2=[4 1 0 2];

nodal_values = 6 ...
                -3*(lege_eval(X(1,:),ll1(1),a,b).*...
                    herm_eval(X(2,:),ll1(2),mu,sig).*...
                    cheb_eval(X(3,:),ll1(3),a2,b2).*...
                    lagu_eval(X(4,:),ll1(4),lambda))...
                +1*(lege_eval(X(1,:),ll2(1),a,b).*...
                    herm_eval(X(2,:),ll2(2),mu,sig).*...
                    cheb_eval(X(3,:),ll2(3),a2,b2).*...
                    lagu_eval(X(4,:),ll2(4),lambda));


domain={[a;b], [mu;sig], [a2;b2], lambda};

[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,{'legendre','hermite','chebyshev','laguerre','generalized laguerre','jacobi'});

sel = abs(modal_coeffs)>1e-12;
[K(sel,:),modal_coeffs(sel,:)]


%% a vector-valued case: each component of the vector-valued field is considered individually

clear

flag='hermite';

% means
mu=[1 -1];
% standard deviations
sig=[2 3];



% The sparse grid. Note that here I don't really care what points I have used to build the nodal
% interpolant. They could even be evenly spaced points! It's just a bunch of evaluations of a function.
% of course bad points make the matrix ill conditioned
N=2; w=5; 
knots1=@(n) knots_normal(n,mu(1),sig(1)); 
knots2=@(n) knots_normal(n,mu(2),sig(2)); 
lev2knots=@lev2knots_lin; 
idxset=@(i) sum(i-1);
S=create_sparse_grid(N,w,{knots1,knots2},lev2knots,idxset);
Sr=reduce_sparse_grid(S);



% Same procedure as before, now with a linear comb of Hermite polynomials
% NODAL_VALUES is a matrix VxM, where V is the output dimensionality and M is the number of points.
% in other words, evaluations of the function-values are stored as columns (same as EVALUATE_ON_SPARSE_GRID)
% so, we no longer need to transpose in nodal_values1 and 2

X=Sr.knots;

nodal_values1 = 6-4*herm_eval_multidim(X,[3 1],mu,sig)+...
    7*herm_eval_multidim(X,[1 2],mu,sig);

nodal_values2 = 0.5-2*herm_eval_multidim(X,[1 1],mu,sig)+...
    5*herm_eval_multidim(X,[2 3],mu,sig);

nodal_values = [nodal_values1; nodal_values2];

domain=[mu;sig];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,flag);

[K,modal_coeffs]
