%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2014 L. Tamellini, F. Nobile
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

S=smolyak_grid(N,w,knots,lev2knots,idxset);
Sr=reduce_sparse_grid(S);

% the domain of the grid
domain=[-ones(1,N); ones(1,N)];


% compute a legendre polynomial over the sparse grid
X=Sr.knots;
nodal_values = 4*lege_eval_multidim(X,[4 0],-1,1)'+...
    2*lege_eval_multidim(X,[1 1],-1,1)';

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
N=2; w=5; knots=@(n) knots_uniform(n,-1,1,'nonprob'); lev2knots=@lev2knots_lin; idxset=@(i) sum(i-1);
[S,C]=smolyak_grid(N,w,knots,lev2knots,idxset,get_interval_map(a,b,'uniform'));
Sr=reduce_sparse_grid(S);



% evaluate a linear comb of Legendre polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*lege_eval_multidim(X,[3 0],a,b)'+...
    7*lege_eval_multidim(X,[1 1],a,b)';


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
knots=@(n) knots_gaussian(n,0,1); 
lev2knots=@lev2knots_lin; 
idxset=@(i) sum(i-1);
S=smolyak_grid(N,w,knots,lev2knots,idxset,get_interval_map(mu,sig,'gaussian'));
Sr=reduce_sparse_grid(S);




% Same procedure as before, now with a linear comb of Hermite polynomials

X=Sr.knots;

nodal_values = 6-4*herm_eval_multidim(X,[3 1],mu,sig)'+...
    7*herm_eval_multidim(X,[1 2],mu,sig)';


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
N=2; w=5; knots=@(n) knots_uniform(n,-1,1,'nonprob'); lev2knots=@lev2knots_lin; idxset=@(i) sum(i-1);
[S,C]=smolyak_grid(N,w,knots,lev2knots,idxset,get_interval_map(a,b,'uniform'));
Sr=reduce_sparse_grid(S);



% evaluate a linear comb of Chebyshev polynomials on the grid: it should be recovered exactly by the
% procedure (with a sufficient number of points in the grid)
X=Sr.knots;

nodal_values = 6-3*cheb_eval_multidim(X,[3 2],a,b)'+...
    7*cheb_eval_multidim(X,[1 1],a,b)';


domain=[a;b];
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'chebyshev');

[K,modal_coeffs]
