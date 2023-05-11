% This tutorial is a "hands-on" manual of the Sparse Grid Matlab Kit a contains examples of use of the main functions.
% More examples can be found in the files TEST_*.m in this folder. A dedicated tutorial for the function 
% ADAPT_SPARSE_GRID can be found in TUTORIAL_ADAPTIVE.m
%
%
%
% This file is structured as follows. Each part is composed by one or more "matlab sections" that can be 
% executed individually by the command "Run Section" (CTRL+ENTER)
% 
%
%
% PART 0: addtopath / set verbosity
%
%
% PART 1: introduction
%   - what is a sparse grid
%   - ingredients of a sparse grid. 1d knots
%   - ingredients of a sparse grid. lev2knots function
%   - ingredients of a sparse grid. multi-index set
%   - data-strucure
%   - modify the domain of a sparse grid
%   - reduce a sparse grid
%
%
% PART 2: evaluate a function on a sparse grid 
%   - basics
%   - use recycling featur
%   - recycle from a "list of points"
%   - use recycling feature for vector output
%   - use parallel feature
%
%
% PART 3: integration on sparse grids - basics
%   - integration on tensor grids 
%   - use other quadrature knots
%   - modify quadrature domain
%   - compute moments of random variables
%   - recycle evaluations from previously computed grids and parallel computation 
%   - how to build more complex sparse grids. anisotropic grids 
%   - how to build more complex sparse grids. use multiidx_box_set
%
%
% PART 4: interpolation on a sparse grid - basics
%   - interpolation on a tensor grid
%   - interpolation error on sparse grid points
%   - plot sparse grid interpolant
%
%
% PART 5: compute the g-pce of a function given its sparse grid approximation
%
%
% PART 6: sparse-grids based sensitivity analysis
%   -  compute Sobol Indices of f
%   -  compute gradients of sparse grid interpolant of f
%
%
% PART 7: export sparse grid on file




%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------




%% PART 0: INSTALL / ADD TO PATH / VERBOSITY

clc
clear
addpath(genpath(pwd)) % do not use addpath(genpath(./)), it won't work properly
disp('path set')


% to suppress most of text output on screen use
%
% global MATLAB_SPARSE_KIT_VERBOSE;
% MATLAB_SPARSE_KIT_VERBOSE = 0;
%
% to resume text output, set
%
% MATLAB_SPARSE_KIT_VERBOSE = 1;


%% PART 1: INTRODUCTION - WHAT IS A SPARSE GRID

% A sparse grid is a linear combination of many tensor grids on R^N (parameter space). 
% Each of the tensor grids included has ``few points''. With suitable linear combinations
% of such grids, it is possible to achieve good accuracy in quadrature and interpolation, 
% with a computational cost lower than using a single tensor grid 

% run these commands to build a sparse grid and visualize each component

N=2; % approximation of two variables
knots=@(n) knots_CC(n,-1,1,'nonprob'); % knots
w = 3; %level
S = create_sparse_grid(N,w,knots,@lev2knots_doubling); % grid

% visualization

% plot the grid itself
plot_sparse_grid(S,[],'color','k','marker','o','MarkerFaceColor','k');

% each component
figure
s_max=length(S);
k=0;
for s=1:s_max
    if ~isempty(S(s).size) % some grids are not included in the linear combination
        k=k+1;
        subplot(2,4,k)
        % we use again plot_sparse grids, which can plot tensor grids too
        plot_sparse_grid(S(s),[],'color','k','marker','o','MarkerFaceColor','k');
        axis square
        %pause
    end
end



%% PART 1: INTRODUCTION - INGREDIENTS OF A SPARSE GRID. 1D KNOTS 

% each of the tensor grids in the sparse grid is built by taking cartesian products of 1D distribution of
% points (in general a different number of points in each direction). The Sparse Grid Matlab Kit provides
% several knots families. These functions also return the quadrature weights associated to the knots
% (more on this later)

% Gauss-Legendre points: quadrature points to approximate integrals like \int_a^b f(x) 1/(b-a) dx with n points
n=5; a=1; b=4;
x=knots_uniform(n,a,b);

figure
plot(x,0*x,'ok','MarkerFaceColor','k','DisplayName','5 GL points')
grid on

% Clenshaw-Curtis points: nested quadrature points to approximate integrals like \int_a^b f(x) 1/(b-a) dx with n
% points. If one "doubles" the number of points, the new points will include the old ones

hold on

n=5; a=1; b=4;
x=knots_CC(n,a,b);
plot(x,1 + 0*x,'or','MarkerFaceColor','r','DisplayName','5 CC points')


n=9; a=1; b=4;
x=knots_uniform(n,a,b);
plot(x,-1 + 0*x,'ob','MarkerFaceColor','b','DisplayName','9 GL points (does NOT includes the 5 points)')

n=9; a=1; b=4;
x=knots_CC(n,a,b);
plot(x,2 + 0*x,'og','MarkerFaceColor','g','DisplayName','9 CC points (includes the 5 points)')


ylim([-1.5 4])
legend show



% Leja points: nested quadrature points to approximate integrals like \int_a^b f(x) 1/(b-a) dx with n
% points. Three different kind of Leja points are available: Line Leja, sym-Line Leja, p-disk Leja (see
% leja_points.m for more details). All Leja points are nested by construction

figure

% ------- line leja ----------
n=5; a=1; b=4;
x=knots_leja(n,a,b,'line');
plot(x,1 + 0*x,'or','MarkerFaceColor','r','DisplayName','5 Line Leja points')

hold on 

n=9; a=1; b=4;
x=knots_leja(n,a,b,'line');
plot(x,2 + 0*x,'or','MarkerFaceColor','r','DisplayName','9 Line Leja points')


% ------- sym leja ----------
n=5; a=1; b=4;
x=knots_leja(n,a,b,'sym_line');
plot(x,3 + 0*x,'ok','MarkerFaceColor','k','DisplayName','5 sym-line Leja points')

hold on 

n=9; a=1; b=4;
x=knots_leja(n,a,b,'sym_line');
plot(x,4 + 0*x,'ok','MarkerFaceColor','k','DisplayName','9 sym-Line Leja points')



% ------- p-disk leja ----------
n=5; a=1; b=4;
x=knots_leja(n,a,b,'p_disk');
plot(x,5 + 0*x,'ob','MarkerFaceColor','b','DisplayName','5 p-Disk Leja points')

hold on 

n=9; a=1; b=4;
x=knots_leja(n,a,b,'p_disk');
plot(x,6 + 0*x,'ob','MarkerFaceColor','b','DisplayName','9 p-Disk Leja points')

grid on
ylim([-1.5 12])
legend show


%% we also provide equispaced (trapezoidal quadrature rule) and midpoint, 
% which should of course be used sparingly as they are not suitable for high-order 
% interpolation (runge phenomenon)

figure

n=5; a=1; b=4;
x=knots_midpoint(n,a,b);
plot(x,1 + 0*x,'or','MarkerFaceColor','r','DisplayName','5 midpoints')

hold on 

n=6; a=1; b=4;
x=knots_trap(n,a,b);
plot(x,1 + 0*x,'ob','MarkerFaceColor','b','DisplayName','6 trapezoidal points')


n=10; a=1; b=4;
x=knots_midpoint(n,a,b);
plot(x,2 + 0*x,'or','MarkerFaceColor','r','DisplayName','10 midpoints')

n=11; a=1; b=4;
x=knots_trap(n,a,b);
plot(x,2 + 0*x,'ob','MarkerFaceColor','b','DisplayName','11 trapezoidal points (nested with 6, 6+5=11)')

n=30; a=1; b=4;
x=knots_midpoint(n,a,b);
plot(x,3 + 0*x,'or','MarkerFaceColor','r','DisplayName','30 midpoints (nested with 10 midpoints)')

grid on
ylim([-0.5 4.5])
legend show


%% quadrature points to approximate integrals with normal pdf, like 
%
% 1/sqrt(2 sig^2 pi) \int_R f(x) e^{ -(x-mi)^2 / (2 sig^2) } dx 


% Gauss-Hermite points 
n=9; mu=0; sig=1;
x=knots_normal(n,mu,sig);

figure
plot(x,0*x,'ok','MarkerFaceColor','k','DisplayName','9 GH points')
grid on

% Genz-Keister / Kronrod - Patterson Nodes : nested quadrature points to approximate integrals as the previous

hold on
n=3; 
x=knots_GK(n,0,1);
plot(x,1 + 0*x,'or','MarkerFaceColor','r','DisplayName','3 GK points')

n=9; 
x=knots_GK(n,0,1);
plot(x, 2 + 0*x,'ob','MarkerFaceColor','b','DisplayName','9 GK points')

% Normal-Leja : nested quadrature points to approximate integrals as the previous

hold on
n=3; 
x=knots_normal_leja(n,mu,sig,'line'); % another option here is 'sym_line'
plot(x,3 + 0*x,'xr','LineWidth',2,'MarkerFaceColor','r','DisplayName','3 Normal-Leja points')

n=9; 
x=knots_normal_leja(n, mu,sig,'line');
plot(x, 4 + 0*x,'xb','LineWidth',2,'MarkerFaceColor','b','DisplayName','9 Normal-Leja points')

ylim([-1.5 7])
legend show


%% quadrature points to approximate integrals with exponential pdf, like 
%
%  \int_[0,+inf] f(x) lambda e^{ -lambda*x } dx 


% Gauss-Laguerre points
n=12; 
lambda = 1; 

x=knots_exponential(n,lambda);

figure
plot(x,0*x,'ok','MarkerFaceColor','k','DisplayName','12 Gauss-Laguerre points')

% and their Leja Counter part

hold on 
x=knots_exponential_leja(n,lambda);
plot(x,1 + 0*x,'or','LineWidth',2,'MarkerFaceColor','r','DisplayName','12 Exponential-Leja points')

grid on
ylim([-0.5 1.5])
legend show

%% quadrature points to approximate integrals with gamma pdf, like 
%
%  \int_[0,+inf] f(x) beta^(alpha+1)/Gamma(alpha+1)*x^alpha*exp(-beta*x) dx 


%  Gauss-Generalized Laguerre points:  with n points
n=12; 
alpha = 1; beta =2;

x=knots_gamma(n,alpha,beta);

figure
plot(x,0*x,'ok','MarkerFaceColor','k','DisplayName','12 Generalized Gauss-Laguerre points')

% and their Leja Counter part

hold on 
x=knots_gamma_leja(n,alpha,beta);
plot(x,1 + 0*x,'or','LineWidth',2,'MarkerFaceColor','r','DisplayName','12 Gamma-Leja points')

grid on
ylim([-0.5 1.5])
legend show


%% quadrature points to approximate integrals with beta pdf, like 
%
%  \int_[x_a, x_b] f(x) Gamma(alpha+beta+2)/(Gamma(alpha+1)*Gamma(beta+1)*(x_b-x_a)^(alpha+beta+1))*(x-x_a)^alpha*(x_b-x)^beta dx 


%  Gauss-Jacobi points: with n points
n=12; 
x_a = 1; x_b = 3;
alpha = -0.5; beta = 0.5; 

x=knots_beta(n,alpha,beta,x_a,x_b);

figure
plot(x,0*x,'ok','MarkerFaceColor','k','DisplayName','12 Gauss-Jacobi points')

% and their Leja Counter part

hold on
n=12; 
x=knots_beta_leja(n,alpha,beta,x_a,x_b,'line'); 
plot(x,1 + 0*x,'or','LineWidth',2,'MarkerFaceColor','r','DisplayName','12 Beta-Leja points')

n=12; 
x=knots_beta_leja(n,alpha,beta,x_a,x_b,'sym_line'); % recommended only if alpha = beta 
plot(x,2 + 0*x,'ob','LineWidth',2,'MarkerFaceColor','b','DisplayName','12 Sym. Beta-Leja points')

grid on
ylim([-0.5 3])
legend show

%% quadrature points to approximate integrals with triangular pdf
% (i.e. a linear decreasing pdf over the interval [a,b]), like 
%
%  \int_[a,b] f(x) * 2/(b-a)^2 * (b-x) dx

n = 12; 
a = 0; b = 2;
x=knots_triangular_leja(n,a,b);

figure
plot(x,0*x,'ok','MarkerFaceColor','k','DisplayName','12 Triangular-Leja points')
grid on
legend show


%% PART 1: INTRODUCTION - INGREDIENTS OF A SPARSE GRID. LEV2KNOTS FUNCTION.

% in view of building sparse grids, it is useful to order quadrature/interpolation rules in sequences, i.e. to
% introduce levels for the rules. The Sparse Grid Matlab Kit provides 5 functions to this end:

% -> lev2knots_lin     
%
% adds 1 point from one level to the next: consecutive quadrature/interpolation rules have
% 1,2,3,4,..points

clc
lev2knots_lin([1 2 3 4 5])


% -> lev2knots_2step     
%
% adds 2 points from one level to the next: consecutive quadrature/interpolation rules have
% 1,3,5,7,..points

lev2knots_2step([1 2 3 4 5])


% -> lev2knots_doubling 
%
% "doubles" the number of points from one level to the next: consecutive rules have 1,3,5,9,17... points

lev2knots_doubling([1 2 3 4 5])


% -> lev2knots_tripling 
%
% triples the number of points from one level to the next: consecutive rules have 1,3,9,27... points.


lev2knots_tripling([1 2 3 4 5])



% -> lev2knots_GK
%
% needed when using GK knots which are tabulated. consecutive rules have 1,3,9,19,35 points. The latter
% is the finest resolution possible

lev2knots_GK([1 2 3 4 5])



%% PART 1: INTRODUCTION - INGREDIENTS OF A SPARSE GRID. MULTI-INDEX SET

% the last ingredient to specify when building a sparse grid is the set of tensor grids to be used. The
% algorithm will then take care of computing the coefficients of the linear combination of these grids
% (note that such coefficients may be 0 as well). 


% The most convenient way to specify tensor grids is to use multi-index notation: every grid is
% associated to a multiindex, that is a vector of integer numbers. 
% Each number in the vector tells the level of the quadrature rule used in each direction
% of the parameter space. E.g. : the multiindex [3 5] is associated to the tensor grid
% using a quad rule of level 3 in direction 1, and level 5 in direction 2. The actual number of points in
% each direction depends by the level-knots relation specified by the lev2knots_*** function.

clc
N=2; 
ii=[3 5];
knots=@(n) knots_uniform(n,-1,1,'nonprob'); % knots

S_lin=tensor_grid(N,lev2knots_lin(ii),knots);
S_doub=tensor_grid(N,lev2knots_doubling(ii),knots);

figure

plot_sparse_grid(S_doub,[],'color','r','marker','s','MarkerFaceColor','r','DisplayName','lev2knots-nested');
hold on
plot_sparse_grid(S_lin,[],'color','k','marker','o','MarkerFaceColor','k','DisplayName','lev2knots-lin');

legend show
set(legend,'Location','SouthOutside')


% there are two ways of specifying the set of multindices to be used. 
%
% 1) The first one is to use the parameters "level" and "idxset" of the function CREATE_sPARSE_GRID. 
% In this case, the multiindex set will include all the multiindices that satisfy the inequality
%
% idxset(ii)<= level
%
% by default, idxset is set to @(ii) sum(ii-1). The combination of idxset function and lev2knots function
% defines the sparse grid type: using @(ii) sum(ii-1) with lev2knots_lin results in the so-called TD
% (Total Degree) tensor grid, while  @(ii) sum(ii-1) with lev2knots_doubling in the original SM (Smolyak) grid.
% Some choices are available by using the function 
%
%  [lev2nodes,idxset] = DEFINE_FUNCTIONS_FOR_RULE(rule,rates)
%
% but any other set satisfying the so-called ``admissibility condition'' 
% (see e.g. Gerstner-Griebel ``Dimension-Adaptive Tensor-Product Quadrature'') can be used.

clc
N=2; 
knots=@(n) knots_uniform(n,-1,1,'nonprob'); 
w = 5; %level

[lev2knots,idxset]=define_functions_for_rule('TD',N);
S_TD = create_sparse_grid(N,w,knots,lev2knots,idxset); % grid

[lev2knots,idxset]=define_functions_for_rule('HC',N);
S_HC = create_sparse_grid(N,w,knots,lev2knots,idxset); % grid

% plot the grid itself
figure
plot_sparse_grid(S_TD,[],'color','k','marker','o','MarkerFaceColor','k');
legend('TD-grid')


figure
plot_sparse_grid(S_HC,[],'color','k','marker','o','MarkerFaceColor','k');
legend('HC-grid')

% 2) The second one is to use the function CREATE_SPARSE_GRID_MULTIIDX_SET, where one specifies exactly
% the set of multiindex that one wishes to use. Again, the set has to satisfy
% the ``admissibility condition'', and the rows have to be in lexicographic order. 

C=[
    1 1;
    1 2;
    1 3;
    1 4;
    2 1;
    2 2;
];

[adm,C] = check_set_admissibility(C);

S_M = create_sparse_grid_multiidx_set(C,knots,lev2knots);

figure
plot_sparse_grid(S_M,[],'color','b','marker','d','MarkerFaceColor','b');
axis([-1 1 -1 1])


%% the package provides two functions to generate multi-index sets. 

% a) MULTIIDX_BOX_SET generates all multiindices jj that are component-wise less than or
% equal to some other index ii. The minimal value of the components of the indices to be generated can be either 0 or 1. For instance

jj=[2 3]; 
C=multiidx_box_set(jj,0);
D=multiidx_box_set(jj,1);

% the package comes with a convenience function to plot a multiidx set

figure
plot_multiidx_set(C,'xr','MarkerFaceColor','r','LineWidth',2,'MarkerSize',12,'DisplayName','Multiidx box set, min=0')
hold on
plot_multiidx_set(D,'ok','MarkerFaceColor','k','DisplayName','Multiidx box set, min=1')
axis([-0.5 4 -0.5 4])
legend show


% b) MULTIIDX_GEN generates the set of all indices ii such that rule(ii)<=w, where rule is a function that takes as input a row vector
% (or a matrix where each multiidx is stored as a row) and returns a scalar value (or a column vector with the result of the operation applied
% to each row of the input index vector). Again, the minimum index can be 0 or 1:

N=2;
w=7;
rule=@(I) sum(I,2); 
E=multiidx_gen(N,rule,w,0);
F=multiidx_gen(N,rule,w,1);


figure 
plot_multiidx_set(E,'xr','MarkerFaceColor','r','LineWidth',2,'MarkerSize',12,'DisplayName','Multiidx gen, min=0')
hold on
plot_multiidx_set(F,'ok','MarkerFaceColor','k','DisplayName','Multiidx gen, min=1')
legend show
axis([-0.5 8 -0.5 8])

% incidentally, PLOT_MULTIIDX_SET works also for N=3. For larger dimensions, one needs to input the subset of dimensions that are to be plotted
G=multiidx_box_set([2 3 5],1);
figure
plot_multiidx_set(G)


H=multiidx_box_set([2 3 5 4],1);
figure
plot_multiidx_set(G(:,[1 3]))


%% when building a large sparse grid, it might be useful to recycle from previous grids to speed-up the computation

clc
clear

knots=@(n) knots_normal(n,0,1);
lev2knots=@lev2knots_lin;

N=20;
w=4;
S=create_sparse_grid(N,w,knots,lev2knots,@(i) prod(i));
w=5;
disp('build grid without recycling')
tic
T=create_sparse_grid(N,w,knots,lev2knots,@(i) prod(i));
toc
tic
T_rec=create_sparse_grid(N,w,knots,lev2knots,@(i) prod(i),S);
toc

isequal(T,T_rec) % sometimes fields like knots or weights might differ at machine precision

%% note that the following call is also valid: 
% T_rec=create_sparse_grid(N,w,knots,lev2knots,@(i) prod(i),[]);
% this is useful in iterative loops like:
clc
tic
for w=1:7
    % build grid
    T=create_sparse_grid(N,w,knots,lev2knots,@(i) prod(i));
    % then do something ...
end
toc


tic
T_old=[];
for w=1:7
    % build grid
    T=create_sparse_grid(N,w,knots,lev2knots,@(i) prod(i),T_old);
    T_old = T;
    % then do something ...
end
toc

%% the same functionality is also available for create_sparse_grid_multiidx_set

clear
clc

knots=@(n) knots_normal(n,0,1);
lev2knots=@lev2knots_lin;
ibox= [3 4 2 4 2];
[~,C] = multiidx_box_set(ibox,1);
D = sortrows([C; 2 5 2 2 6]);

S=create_sparse_grid_multiidx_set(C,knots,lev2knots);

tic
T=create_sparse_grid_multiidx_set(D,knots,lev2knots);
toc
tic
T_rec = create_sparse_grid_multiidx_set(D,knots,lev2knots,S);
toc
isequal(T,T_rec)



%% PART 1: INTRODUCTION - DATA-STRUCURE

% A sparse grid is represented as a vector of structures. Each element is a tensor grid, with fields
% containing the knots, the corresponding integration weights, its coefficient in the linear combination,
% and the number of points.

% In general, the following conventions hold:
% -> points in the space of parameters are columns-vector
% -> multiindices are row-vector


%% PART 1: INTRODUCTION - MODIFY THE DOMAIN OF A SPARSE GRID

% it is easy to modify the domain of a sparse grid from (-1,1)^N to other hyper-rectangles.

clc
clear
N=2;

% generate knots on the desired hyper-rectangle (here (0,2)^2 )
knots=@(n) knots_CC(n,0,2,'nonprob');
w = 4;
S = create_sparse_grid(N,w,knots,@lev2knots_doubling);



figure
plot_sparse_grid(S);
legend('grid S')
set(legend,'Location','NorthEastOutside')

% one can mix different intervals / different knots families on different directions. 

clc
clear
N=2;

knots1=@(n) knots_CC(n,0,2,'nonprob');
knots2=@(n) knots_uniform(n,-1,5,'nonprob');
w = 4;
S = create_sparse_grid(N,w,{knots1,knots2},{@lev2knots_doubling,@lev2knots_lin});

figure
plot_sparse_grid(S);



%% PART 1: INTRODUCTION - REDUCE A SPARSE GRID

% The tensor grids forming the sparse grid may have points in common (even when using non-nested points).
% To save computational time during e.g. evaluation of a function on a sparse grid, it is then important
% to get rid of these repetions. To this end, use the function reduce_sparse_grid. The quadrature weights
% are of course consistently modified. The field "size" tells the number in the reduced grid

clc
clear 
N=2;
w=5;
knots=@(n) knots_CC(n,-1,1,'nonprob'); 


[lev2nodes,idxset] = define_functions_for_rule('SM',N); 
S = create_sparse_grid(N,w,knots,lev2nodes,idxset); 
Sr=reduce_sparse_grid(S);


fprintf('size of original grid: %i\n',size([S.knots],2))
fprintf('size of reduced  grid: %i\n',size(Sr.knots,2))
fprintf('Sr.size: %i\n',Sr.size)


figure
subplot(1,2,1)
plot_sparse_grid(S,[],'color','b','marker','d','MarkerFaceColor','b');
axis square
legend('original grid')
set(legend,'Location','SouthOutside')

subplot(1,2,2)
plot_sparse_grid(Sr,[],'color','b','marker','d','MarkerFaceColor','b');
axis square
legend('reduced grid')
set(legend,'Location','SouthOutside')

%% The Kit provides a short-hand to create and reduce a "vanilla sparse grid", i.e. 
%  - Clenshaw--Curtis points in [-1,1] 
%  - lev2knots_doubling 
%  - multi-index set: sum(ii-1) \leq w 
%
%  (cf define_functions_for_rule('SM') )


clear 
clc
N = 2;
w = 3;
[S,Sr] = create_sparse_grid_quick_preset(N,w);


%% PART 2: EVALUATE A FUNCTION ON A SPARSE GRID - BASICS

% the kit comes with the function evaluate_on_sparse_grid, that allows to  evaluate a function on the points of a sparse grid, and provides
% -> recycling of previous evaluations 
% -> support for parallel evaluations. 
% Works for scalar-valued as well as vector-valued functions. Sparse grids passed as input must be reduced


clc
clear

fs=@(x) sum(x);
fv=@(x) 2*x;

N=2; w=3;
S=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);


% plain use of evaluate_on_sparse_grid: no recycling, no parallel
evals_plain_fs=evaluate_on_sparse_grid(fs,Sr);
evals_plain_fv=evaluate_on_sparse_grid(fv,Sr);

% a direct computation
pts = Sr.size;

os=size(fs(Sr.knots(:,1)),1);
ov=size(fv(Sr.knots(:,1)),1);

evals_direct_fs = zeros(os,pts);
evals_direct_fv = zeros(ov,pts);

for i=1:pts
    evals_direct_fs(:,i)=fs(Sr.knots(:,i));
    evals_direct_fv(:,i)=fv(Sr.knots(:,i));
end



% compare the two values
find(evals_plain_fs~=evals_direct_fs)
find(evals_plain_fv~=evals_direct_fv)

%% PART 2: EVALUATE A FUNCTION ON A SPARSE GRID - USE RECYCLING FEATURE

clc
clear

f=@(x) sum(x);

N=2; w=3;
S=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=4;
T=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);

evals_non_rec=evaluate_on_sparse_grid(f,Tr);
evals_rec=evaluate_on_sparse_grid(f,T,Tr,evaluate_on_sparse_grid(f,Sr),S,Sr);

max(abs(evals_non_rec(:)-evals_rec(:))) 

%% PART 2: EVALUATE A FUNCTION ON A SPARSE GRID - RECYCLE FROM A "LIST OF POINTS"

% it is also possible to recycle from a list of points. However, the algorithm used to detect
% which points are to be evaluated is much slower than the previous case for N large

clc
clear

f=@(x) sum(x);

N=20; w=1;
S=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=2;
T=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);

evals_non_rec=evaluate_on_sparse_grid(f,Tr);
tic
evals_rec=evaluate_on_sparse_grid(f,T,Tr,evaluate_on_sparse_grid(f,Sr),S,Sr);
toc
% pretend we only know the list of points Sr.knots, to see the difference in performance ... 
tic
evals_rec_slow=evaluate_on_sparse_grid(f,T,Tr,evaluate_on_sparse_grid(f,Sr),[],Sr.knots);
toc
max(abs(evals_non_rec(:)-evals_rec(:))) 
isequal(evals_rec,evals_rec_slow)




%% PART 2: EVALUATE A FUNCTION ON A SPARSE GRID - USE RECYCLING FEATURE FOR VECTOR OUTPUT

clc
clear

f=@(x) 2*x;

N=2; w=1;
S=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=2;
T=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);

evals_non_rec=evaluate_on_sparse_grid(f,Tr);
evals_rec=evaluate_on_sparse_grid(f,T,Tr,evaluate_on_sparse_grid(f,Sr),S,Sr);

max(abs(evals_non_rec(:)-evals_rec(:))) 

%% PART 2: EVALUATE A FUNCTION ON A SPARSE GRID - USE PARALLEL FEATURE

% parallel computation can be used both with and without recycling. The parallel procedure gets activated only
% if at least X evaluations are queried,  with X specified by the user. This is because parallel computations have some
% communication overhead, therefore if function evaluations are fast the parallel evaluation may actually result slower
% than the serial.


clc
clear

f=@(x) sum(x);

N=2; w=3;
S=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=4;
T=create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);

if ~check_if_parallel_on()
    activate_parallel() % optional argument to specify how many workers
end
X=0;
evals_1=evaluate_on_sparse_grid(f,T,Tr,[],[],[],X);
X=10;
evals_2=evaluate_on_sparse_grid(f,T,Tr,[],[],[],X);
X=100;
evals_3=evaluate_on_sparse_grid(f,T,Tr,evaluate_on_sparse_grid(f,Sr),S,Sr,X,1e-14);

figure
plot(evals_1)
hold on
plot(evals_2,'x')
plot(evals_3,'o')



if check_if_parallel_on()
    close_parallel()
end

% this command will now throw an error
%
%   evals_1=evaluate_on_sparse_grid(f,Tr,[],[],0);
%
% Error using evaluate_on_sparse_grid>simple_evaluate (line 179)
% no open matlabpool session detected


%% PART 2: EVALUATE A FUNCTION ON A SPARSE GRID
%
% see test_evaluate_on_sparse_grids.m for more examples 


%% PART 3: INTEGRATION ON SPARSE GRIDS - BASICS

% In this part we show how to use the Kit to perform high-dimensional quadrature. We consider the
% following function, for which we know the analytic expression of the integral
%
%   f(x) = prod(1/sqrt(x_i + b))  in [-1,1]^N

clc
clear
f = @(x,b) prod(1./sqrt(x+b));
b=3;
N = 4;
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b));
I_ex = I_1d^N;

% generate the knots and the SM grid. 'nonprob' means we are integrating w.r.t. the pdf rho(x)=1 and not rho(x)=1/prod(b_i - a_i)
knots=@(n) knots_CC(n,-1,1,'nonprob');
w = 4;
S = create_sparse_grid(N,w,knots,@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

% compute integral
I=f([Sr.knots],b)*[Sr.weights]'  %#ok<NOPTS>

% alternatively use
I2=quadrature_on_sparse_grid(@(x)f(x,b) , Sr); % Sr must be reduced here

disp('----------')
disp('difference between values')

I-I2  %#ok<MNEFF,NOPTS>

% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)

% sometimes, we have access to the evaluations of f from earlier code, then we just need to do the linear combination. The package provides
% a convenience wrapper to this end,  instead of typing f_vals*Sr.weights'  

disp('----------')
disp('convenience wrapper')

f_vals = f([Sr.knots],b);
I3 = quadrature_on_sparse_grid(f_vals,Sr);

I2-I3

% the convenience wrapper can also handle the case of computing quadrature for multiple functions at the same time. The values of each function
% must be stored as rows of a matrix
many_f = [f_vals; f_vals; f_vals; f_vals; f_vals];
quadrature_on_sparse_grid(many_f,Sr)


%% PART 3: INTEGRATION ON TENSOR GRIDS 

% integration on tensor grids 

clc
clear
f = @(x,b) prod(1./sqrt(x+b));
b = 3;
N = 2;
I_1d = (2*sqrt(1+b)-2*sqrt(-1+b));
I_ex = I_1d^N;

% let's build a tensor grid with the following choices
knots=@(n) knots_CC(n,-1,1,'nonprob');
ii = [6,5]; 
lev2knots = @lev2knots_doubling;

% create the tensor grid, convert it to sparse
T = tensor_grid(N,lev2knots(ii),knots);  
S = tensor_to_sparse(T);
% with the call above, S.idx is automatically set to the number of points in each dir, 
% which implies S.idx = lev2knots(ii). If you want S.idx = ii, use the following call
% S = tensor_to_sparse(T,ii);
Sr = reduce_sparse_grid(S);
is_sparse_grid(S)

I=quadrature_on_sparse_grid(@(x)f(x,b) , Sr); % Sr must be reduced here

% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)

%% PART 3: INTEGRATION - USE OTHER QUADRATURE KNOTS

% as already seen in the introduction, other quadrature knots are available

clc
clear

f = @(x,b) prod(1./sqrt(x+b));
b=3;
N = 4;
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b));
I_ex = I_1d^N;


knots=@(n) knots_uniform(n,-1,1,'nonprob');
w = 4;
S = create_sparse_grid(N,w,knots,@lev2knots_doubling);

Sr=reduce_sparse_grid(S);

I=quadrature_on_sparse_grid(@(x)f(x,b) , Sr);

% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)



%% PART 3: INTEGRATION - MODIFY QUADRATURE DOMAIN
clear
clc

% suppose integrating over (-1,3)^N
f = @(x,b) prod(1./sqrt(x+b));
b=3;
N = 4;

I_1d=(2*sqrt(3+b)-2*sqrt(-1+b));
I_ex = I_1d^N;

% generate knots in (-1,3)
knots=@(n) knots_CC(n,-1,3,'nonprob');
w = 6;
S = create_sparse_grid(N,w,knots,@lev2knots_doubling);
Sr= reduce_sparse_grid(S);

I=quadrature_on_sparse_grid(@(x)f(x,b) , Sr);


% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)


%% PART 3: INTEGRATION - COMPUTE MOMENTS OF RANDOM VARIABLES

% here we compute E[f(x)] = \int_{[-2 1]x[0.5 6]} f(x) 1/(3*5.5) dx, (3*5.5 is the size of the domain)
clc
clear

f = @(x,b) prod(1./sqrt(x+b));
b=3;
N = 2;

I_ex = 1/3/5.5*(2*sqrt(1+b)-2*sqrt(-2+b))*(2*sqrt(6+b)-2*sqrt(0.5+b))  %#ok<NOPTS>

% the best-practice is to generate knots on (-2,1) and (0.5,6), specifying 'prob' as input to the
% knots-generatic function
knots1=@(n) knots_CC(n,-2,1,'prob'); % knots1=@(n) knots_CC(n,-2,1); would work as well as 'prob' is the default value
knots2=@(n) knots_CC(n,0.5,6,'prob'); % knots2=@(n) knots_CC(n,0.5,6); would work as well as 'prob' is the default value
w = 6;
S = create_sparse_grid(N,w,{knots1,knots2},@lev2knots_doubling);
Sr = reduce_sparse_grid(S);
I=quadrature_on_sparse_grid(@(x)f(x,b) , Sr);


% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)



%% PART 3: INTEGRATION - RECYCLE EVALUATIONS FROM PREVIOUSLY COMPUTED GRIDS AND PARALLEL COMPUTATION

% just as evaluate_on_sparse_grid, quadrature_on_sparse_grid provides evaluation recycling and
% parallel evaluation
 
clc
clear

f = @(x,b) prod(1./sqrt(x+b));
b=5;
N = 2;

% the starting grid
w = 7;
S = create_sparse_grid(N,w,@(n) knots_CC(n,-2,1,'prob'),@lev2knots_doubling);
Sr = reduce_sparse_grid(S);
[IS,evals_S]=quadrature_on_sparse_grid(@(x)f(x,b), Sr);

% the new grid
w = 8;
T = create_sparse_grid(N,w,@(n) knots_CC(n,-2,1,'prob'),@lev2knots_doubling);
Tr = reduce_sparse_grid(T);
% the recycling call. Other optional arguments can turn on parallel
% evaluation and tune the tolerance for two points to be considered equal.
% See help quadrature_on_sparse_grid and test_evaluate_on_sparse_grid
IT_rec=quadrature_on_sparse_grid(@(x)f(x,b),T,Tr,evals_S,S,Sr);


% the non-recycling call
IT=quadrature_on_sparse_grid(@(x)f(x,b) , Tr);

if ~check_if_parallel_on()
    activate_parallel() % optional argument to specify how many workers
end

% the parallel call with no recycling
IT2= quadrature_on_sparse_grid(@(x)f(x,b) , T, Tr, [],[],[],0);

% the parallel call with recycling
IT3=quadrature_on_sparse_grid(@(x)f(x,b),T,Tr,evals_S,S,Sr,0);


disp('-------------')
disp('difference between the results')
[IT_rec; IT; IT2; IT3] %#ok<NOPTS>

%% PART 3: INTEGRATION - HOW TO BUILD MORE COMPLEX SPARSE GRIDS. ANISOTROPIC GRIDS 

clear
clc

f = @(x,b) prod(1./sqrt(x+b)); 
b=3; 
N = 4; 
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b)); 
I_ex = I_1d^N;


% specify a rule like in Back Nobile Tamellini Tempone, `Stochastic Spectral Galerkin and Collocation...a  numerical comparison''
rates=[1 2 2 2];
knots=@(n) knots_uniform(n,-1,1,'nonprob');
[lev2nodes,idxset] = define_functions_for_rule('TD',rates); 
w=4;
[S2] = create_sparse_grid(N,w,knots,lev2nodes,idxset); 

 
% use it to compute integral 
I=quadrature_on_sparse_grid(@(x) f(x,b),reduce_sparse_grid(S2));

% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)


%% PART 3: INTEGRATION - HOW TO BUILD MORE COMPLEX SPARSE GRIDS. USE MULTIIIDX_BOX_SET

% As seen in the introduction, specify directly the set of multiindices involved. 
% Here, we generate the box set of all multiindices <= of [3 5  2 3] in lexicographic order

clc
clear 

f = @(x,b) prod(1./sqrt(x+b)); 
b=3; 
N = 4; 
I_1d=(2*sqrt(1+b)-2*sqrt(-1+b)); 
I_ex = I_1d^N;


C = multiidx_box_set([3 5 2 3],1); % X is C without [3 5 2 3]
knots=@(n) knots_uniform(n,-1,1,'nonprob');

S3=create_sparse_grid_multiidx_set(C,knots,@lev2knots_lin);



% use it to compute integral (-1,1 Lebesgue measure)
I=quadrature_on_sparse_grid(@(x) f(x,b),reduce_sparse_grid(S3));

% compare with exact value
disp('----------')
disp('quad error')
abs(I-I_ex)

%% PART 3: INTEGRATION - CONVERGENCE STUDY

% see test_sparse_quadrature.m



%% PART 4: INTERPOLATION ON A SPARSE GRID - BASICS

% the sparse grid also provides an interpolant / surrogate model for the original function. The
% interpolant can be evaluated in non-grid points. 
%
% All the previous topics (changing the domain, building anisotropic grids ...) apply immediately to 
% the interpolation case. 

clc
clear
f = @(x,b) prod(1./sqrt(x+b)); b=3; N = 4; 

w=8;
knots=@(n) knots_uniform(n,-1,1,'nonprob');
[S] = create_sparse_grid(N,w,knots,@lev2knots_lin); 

Sr=reduce_sparse_grid(S);

%non_grid_points=rand(N,100); 
non_grid_points=[0.5*ones(N,1), zeros(N,1)];

function_on_grid=f(Sr.knots,b);

f_values = interpolate_on_sparse_grid(S,Sr,function_on_grid,non_grid_points);

disp('----------')
disp('Interpolation error')
max( abs( f_values-f(non_grid_points,b) ) )

%% PART 4: INTERPOLATION ON A TENSOR GRID 

% clc
clear
f = @(x,b) prod(1./sqrt(x+b)); 
b=3; N = 2; 

% let's build a tensor grid with the following choices
idx = [10 8];
knots=@(n) knots_uniform(n,-1,1,'nonprob');
lev2knots = @lev2knots_lin;


% create the tensor grid, convert it to sparse
T = tensor_grid(N,lev2knots(idx),knots);  
S = tensor_to_sparse(T);
Sr = reduce_sparse_grid(S);


%non_grid_points=rand(N,100); 
non_grid_points=[0.5*ones(N,1), zeros(N,1)];

function_on_grid=f(Sr.knots,b);

f_values = interpolate_on_sparse_grid(S,Sr,function_on_grid,non_grid_points);

disp('----------')
disp('Interpolation error')
max( abs( f_values-f(non_grid_points,b) ) )


%% PART 4: INTERPOLATION IN 1D

% for 1D interpolation, a dedicated function exists. It's called univariate_interpolant and can operate on
% vector-valued function, like here below, where we interpolate a function with two components, i.e., F: R -> R^2

clear
 
% the two components of f
f1 = @(x) x.^3; 
f2 = @(x) sin(2*x); 

% interpolation points and values
x_interp = linspace(-1,2,4);
F_interp = [f1(x_interp);
            f2(x_interp)];


% evaluate the interpolant on a much finer grid
x_eval = -1:0.01:2;
F_eval_interp = univariate_interpolant(x_interp,F_interp,x_eval);  

% plot
subplot(1,2,1)
plot(x_eval,F_eval_interp(1,:),'DisplayName','interpolant')
hold on
plot(x_interp,F_interp(1,:),'o','DisplayName','interpolation points')
plot(x_eval,f1(x_eval),'DisplayName','true fun')
legend show

subplot(1,2,2)
plot(x_eval,F_eval_interp(2,:),'DisplayName','interpolant')
hold on
plot(x_interp,F_interp(2,:),'o','DisplayName','interpolation points')
plot(x_eval,f2(x_eval),'DisplayName','true fun')
%legend show





%% PART 4: INTERPOLATION ON A SPARSE GRID -  INTERPOLATION ERROR ON SPARSE GRID POINTS

% since the sparse grid is a linear combination of several tensor grid interpolants, the interpolation
% error in a point of the sparse grid is not necessarily zero, unless all tensor interpolats
% include that point


clc
clear
f = @(x,b) prod(1./sqrt(x+b)); b=3; N = 4; 

% a sparse grid with non-nested points: interpolation error in sparse grid points will be
% non-zero in general

w=4;
knots=@(n) knots_uniform(n,-1,1,'nonprob');
S = create_sparse_grid(N,w,knots,@lev2knots_lin); 
Sr=reduce_sparse_grid(S);

non_grid_points=zeros(N,1); 
function_on_grid=evaluate_on_sparse_grid(@(x) f(x,b), Sr);

f_values = interpolate_on_sparse_grid(S,Sr,function_on_grid,non_grid_points);

disp('----------')
disp('Interpolation error - non-nested grid')
max( abs( f_values-f(non_grid_points,b) ) )

% the interpolation error will instead be zero if we use nested points and
% consider e.g. [0 0 0 0] which belongs to all of the tensor grids

knots=@(n) knots_CC(n,-1,1,'nonprob');
T = create_sparse_grid(N,w,knots,@lev2knots_doubling); 
Tr=reduce_sparse_grid(T);

non_grid_points=zeros(N,1); 
function_on_grid=f(Tr.knots,b);

f_values = interpolate_on_sparse_grid(T,Tr,function_on_grid,non_grid_points);

disp('----------')
disp('Interpolation error - nested grid')
max( abs( f_values-f(non_grid_points,b) ) )


%% PART 4: INTERPOLATION ON A SPARSE GRID -  PLOT SPARSE GRIDS INTERPOLANT - case N=2

clear

% define sparse grid over [4,6] x [1,5]
N=2;
aa=[4 1];
bb=[6 5];

% the function to be interpolated
f=@(x) 1./(1+0.5*sum(x.^2)); 


% create a sparse grid and evaluate the function on it
domain = [aa; bb];
knots1=@(n) knots_CC(n,aa(1),bb(1),'nonprob');
knots2=@(n) knots_CC(n,aa(2),bb(2),'nonprob');
w = 4;
S = create_sparse_grid(N,w,{knots1,knots2},@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

values_on_grid=evaluate_on_sparse_grid(f,Sr);

% the plot: several examples of usage 
figure
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid);
view([200 16])

figure
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid,'with_f_values');

figure
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid,'nb_plot_pts',10);

% access to plot handles for further editing is available. E.g., this sets dots to black
figure
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid,'with_f_values','nb_plot_pts',10);
h=gcf;
axes_h = get(h,'Children');
objs_h = get(axes_h,'Children');
set(objs_h(1),'MarkerFaceColor','k');

figure
plot_sparse_grid(Sr)
axis square


%% as expected, the interpolant might be bad if equispaced point are used

clear
clc

f = @(x) 1./(1+(5*x(1)).^2)*1./(1+(5*x(2)).^2); 

a=-1; b=1;
domain = [-1 -1; 1 1];

N= 2; 
w= 6;
lev2knots = @lev2knots_doubling;

% a bad choice: equispaced points, we get them by trap-rule. Build a grid with them
knots_bad = @(n) knots_trap(n,a,b,'nonprob');
S_bad = create_sparse_grid(N,w,knots_bad,lev2knots);
Sr_bad = reduce_sparse_grid(S_bad);
f_values_bad = evaluate_on_sparse_grid(f,Sr_bad);

% a good choice: CC points. Build another grid with them, to compare. Note that the two grids will have the same
% number of points
knots_ok = @(n) knots_CC(n,a,b,'nonprob');
S_ok = create_sparse_grid(N,w,knots_ok,lev2knots);
Sr_ok = reduce_sparse_grid(S_ok);
f_values_ok = evaluate_on_sparse_grid(f,Sr_ok);


% a subplot figure where we compare the two grids and the two interpolants
figure

subplot(2,2,1)
plot_sparse_grid(Sr_bad,[],'o','MarkerSize',4,'LineWidth',2)
axis square

subplot(2,2,2)
plot_sparse_grids_interpolant(S_bad,Sr_bad,domain,f_values_bad,'with_f_values','nb_plot_pts',40);

subplot(2,2,3)
plot_sparse_grid(Sr_ok,[],'o','MarkerSize',4,'LineWidth',2)
axis square

subplot(2,2,4)
plot_sparse_grids_interpolant(S_ok,Sr_ok,domain,f_values_ok,'with_f_values','nb_plot_pts',40);


%% case N=3

clear

% define sparse grid over [4,6] x [1,5] x [2 3]
N=3;
aa=[4 1 2];
bb=[6 5 3];

% the function to be interpolated
f=@(x) 1./(1+0.5*sum(x.^2)); 


% create a sparse grid and evaluate the function on it
domain = [aa; bb];
knots1=@(n) knots_CC(n,aa(1),bb(1),'nonprob');
knots2=@(n) knots_CC(n,aa(2),bb(2),'nonprob');
knots3=@(n) knots_CC(n,aa(3),bb(3),'nonprob');
w = 4;
S = create_sparse_grid(N,w,{knots1,knots2,knots3},@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

values_on_grid=evaluate_on_sparse_grid(f,Sr);

figure
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid);
 
% specify PlotSpec
figure
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid,'with_f_values','nb_plot_pts',10,'nb_contourfs',10,'nb_contourf_lines',40);

% we have two ways of plotting the sparse grid:

% way 1): 3d-plot
figure
plot3_sparse_grid(Sr,[],'o','MarkerSize',8,'MarkerFaceColor','r')
axis square


% way 2): two-dimensional projections (they will look identical in this case)
figure
subplot(1,3,1)
plot_sparse_grid(Sr,[1 2])
axis square
subplot(1,3,2)
plot_sparse_grid(Sr,[2 3])
axis square
subplot(1,3,3)
plot_sparse_grid(Sr,[1 3])
axis square





%% case N>3

clear

% define sparse grid over [-1,1]^7
N=7;
aa=-1*ones(1,N);
bb=ones(1,N);

% the function to be interpolated
f=@(x) 1./(1+0.5*x(1,:).^2+0.25*x(2,:).^2+5*x(3,:).^2+2*x(4,:).^2+0.001*x(5,:).^2+10*x(6,:).^2 + 10*x(7,:).^2); 
%f=@(x) 1./(1+0.5*x(1,:).^2+0.5*x(2,:).^2+0.5*x(3,:).^2+0.5*x(4,:).^2+0.5*x(5,:).^2+0.5*x(6,:).^2 + 0.5*x(7,:).^2); 


% create a sparse grid and evaluate the function on it
domain = [aa; bb];
knots=@(n) knots_CC(n,-1,1,'nonprob');
w = 6;
S = create_sparse_grid(N,w,knots,@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

values_on_grid=evaluate_on_sparse_grid(f,Sr);

% add f_values. Note that there are possibly several points which share the values of the coordinates in the cuts,
% therefore there will be points not on the surface. This helps understanding the fluctuations of the function
% when the coordinates not in the cut are not fixed to their average value. In this specific example, changing the
% values of the frozen variables from their averages happens to lower the value of the function. The function generates
% one new figure per cut
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid,'with_f_values');

% specify cuts. Again, because we are specifying cuts, a new figure per cut is generated.  The code below generates two figures 
% (the first one is empty)
plot_sparse_grids_interpolant(S,Sr,domain,values_on_grid,'two_dim_cuts',[1 4 2 7]);



% we have two ways of plotting the sparse grid:

% way 1): 3d-plot projections, e.g.
figure
plot3_sparse_grid(Sr,[1 2 4],'o','MarkerSize',8,'MarkerFaceColor','r')
axis square


% way 2): two-dimensional projections 
figure
plot_sparse_grid(Sr,[1 2])
axis square




%% PART 4: INTERPOLATION ON A SPARSE GRID - CONVERGENCE STUDY

% see test_sparse_interpolation.m


%% PART 5: COMPUTE THE g-PCE OF A FUNCTION GIVEN ITS SPARSE GRID APPROXIMATION

% the kit provides a function to compute the generalized Polynomial Cahos Expansion (g-PCE) of a function
% of several variables, i.e. the expansion of f in terms of a sum of orthonormal polynomials.
%
% Supported random variables - orthonormal polynomials are
% uniform       | Legendre, Chebyshev
% normal        | Hermite
% exponential   | Laguerre
% gamma         | Generalized Laguerre
% beta          | Jacobi
% 
% The coefficients of these expansions are defined as suitable integrals over the space of parameters, and
% could thus be approximated with sparse grid quadrature. However, a more efficient technique can be
% applied, and it actually implemented in the Kit. It consists in rearranging the sparse grid
% interpolant, which is a linear combination of Lagrange polynomials, as a summation of orthonormal
% polynomials (i.e. performing a change of base to express the same polynomial). Given the relations
% between sparse grids and orthogonal expansion, it is always possible to tune the sparse grid so to
% obtain the gPCE in a precise polynomial space. 
%
% See e.g. Back Nobile Tamellini Tempone, `Stochastic Spectral Galerkin and Collocation...a  numerical
% comparison'' for more details on the sparse grid/orthogonal expansion relation and Tamellini ph.D.
% thesis, chap.6 or MOX report 13/2012 by Formaggia Guadagnini Imperiali Lever Porta Riva Scotti Tamellini
% for details on the conversions
%
% more examples with different kinds of random variables / orthogonal polynomials can be found in test_convert_to_modal.m

clc
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
nodal_values = 4*lege_eval_multidim(X,[4 0],-1,1)+ 2*lege_eval_multidim(X,[1 1],-1,1);

% conversion from the points to the legendre polynomial. I should recover it exactly
[modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,'legendre');

[K,modal_coeffs] %#ok<NOPTS>


%% PART 6: SPARSE-GRIDS-BASED SENSITIVITY ANALYSIS - COMPUTE SOBOL INDICES OF A FUNCTION

clear

% define sparse grid 
aa=[-1 -1 -1];
bb=[ 1  1  1];

% we consider different functions, that we evaluate in one go:
f1=@(x) 1 + x(1,:).^2   + x(2,:).^2 + x(3,:).^2;        
f2=@(x) 1 + 5*x(1,:).^2 + x(2,:).^2 + x(3,:).^2; 
f3=@(x) 1./(1 + x(1,:).^2   + x(2,:).^2 + x(3,:).^2); 
f4=@(x) 1./(1 + 5*x(1,:).^2 + x(2,:).^2 + x(3,:).^2); 

f=@(x) [f1(x); f2(x); f3(x); f4(x)];

% We expect to see these results:
%   f1: has no mixed effects, so the principal and total Sobol indices are identical. Also, it's isotropic, so the indices of each variable are identical
%   f2: no mixed effects as f1, but y_1 contributes more to the variability of f so it has a larger Sobol total/principal index
%   f3: this function has mixed effects (partial derivatives are nonzero),  so the principal and total Sobol index will be different, but equal among random variables
%   f4: mixed effects, and y_1 contributes more to the variability of f so it has larger Sobol indices

% generate a sparse grid
domain = [aa; bb;];
knots=@(n) knots_CC(n,-1,1,'nonprob');
N = length(aa);
w = 5;
S = create_sparse_grid(N,w,knots,@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

values_on_grid=evaluate_on_sparse_grid(f,Sr);

% compute Sobol indices. The function uses internally the function CONVERT_TO_MODAL and it uses the same inputs. It works for scalar-values functions only, 
% so we need to run it 4 times
[Sob_i1,Tot_Sob_i1,Mean1,Var1] = compute_sobol_indices_from_sparse_grid(S,Sr,values_on_grid(1,:),domain,'legendre');
[Sob_i2,Tot_Sob_i2,Mean2,Var2] = compute_sobol_indices_from_sparse_grid(S,Sr,values_on_grid(2,:),domain,'legendre');
[Sob_i3,Tot_Sob_i3,Mean3,Var3] = compute_sobol_indices_from_sparse_grid(S,Sr,values_on_grid(3,:),domain,'legendre');
[Sob_i4,Tot_Sob_i4,Mean4,Var4] = compute_sobol_indices_from_sparse_grid(S,Sr,values_on_grid(4,:),domain,'legendre');

% results are as expected
disp('      f1   |    f2    |   f3    |    f4   ')
disp('Principal Sobol indices')
disp([Sob_i1 Sob_i2 Sob_i3 Sob_i4])
disp('Total Sobol indices')
disp([Tot_Sob_i1 Tot_Sob_i2 Tot_Sob_i3 Tot_Sob_i4])

%% PART 6: SPARSE-GRIDS-BASED SENSITIVITY ANALYSIS - COMPUTE GRADIENTS OF A SPARSE GRID INTERPOLANT (by finite differences)

clear

% define sparse grid over [4,6] x [1,5]
N=2;
aa=[4 1];
bb=[6 5];

% the function to be interpolated and its derivatives
f=@(x) 1./(1+0.5*sum(x.^2)); 
df1 = @(x) -1./((1+0.5*sum(x.^2)).^2)*2*0.5.*x(1,:);
df2 = @(x) -1./((1+0.5*sum(x.^2)).^2)*2*0.5.*x(2,:);



% create a sparse grid and evaluate the function on it
domain = [aa; bb];
knots1=@(n) knots_CC(n,aa(1),bb(1),'nonprob');
knots2=@(n) knots_CC(n,aa(2),bb(2),'nonprob');
w = 4;
S = create_sparse_grid(N,w,{knots1,knots2},@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

values_on_grid=evaluate_on_sparse_grid(f,Sr);

% generate M random points in the domain where we evaluate the derivative of the sparse grid 
% and the true derivative, to check error
M=100;
% use get interval map to go from [-1,1]^N to actual domain
my_map=get_interval_map(aa,bb,'uniform');
eval_points = my_map(rand(N,M)*2-1);


% compute values with function
Grads = derive_sparse_grid(S,Sr,values_on_grid,domain,eval_points);


% error and visualization

max(abs(Grads(1,:) - df1(eval_points)))
max(abs(Grads(2,:) - df2(eval_points)))

figure
hold on; 
plot(Grads(1,:),'-o','DisplayName','Finite Diff'); 
plot(df1(eval_points),'-','DisplayName','true val')
legend show
grid on

figure
hold on; 
plot(Grads(2,:),'-o','DisplayName','Finite Diff'); 
plot(df2(eval_points),'-','DisplayName','true val')
legend show
grid on



%% h is computed automatically in each direction as (b-a)/1E5, but can be adjusted if needed. 
% In the example below, the length of the interval along direction 1 is O(1E-5) so choosing
% the default h would lead to h = O(1E-10), which incurs in numerical cancellations.
% Thus, setting manually a larger value for h helps in reducing the error

N=2;

aa=[4E-5 1];
bb=[6E-5 5];

% the function to be interpolated and its derivatives
f=@(x) 1./(1+0.5*sum(x.^2)); 
df1 = @(x) -1./((1+0.5*sum(x.^2)).^2)*2*0.5.*x(1,:);
df2 = @(x) -1./((1+0.5*sum(x.^2)).^2)*2*0.5.*x(2,:);



% create a sparse grid and evaluate the function on it
domain = [aa; bb];
knots1=@(n) knots_CC(n,aa(1),bb(1),'nonprob');
knots2=@(n) knots_CC(n,aa(2),bb(2),'nonprob');
w = 5;
S = create_sparse_grid(N,w,{knots1,knots2},@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

values_on_grid=evaluate_on_sparse_grid(f,Sr);

% generate M random points in the domain where we evaluate the derivative of the sparse grid 
% and the true derivative, to check error
M=100;
% use get interval map to go from [-1,1]^N to actual domain
my_map=get_interval_map(aa,bb,'uniform');
eval_points = my_map(rand(N,M)*2-1);


% compute values with function
Grads_def = derive_sparse_grid(S,Sr,values_on_grid,domain,eval_points);
h=[1E-7 1E-5];
Grads_man = derive_sparse_grid(S,Sr,values_on_grid,domain,eval_points,h);

% error and visualization

figure
hold on; 
plot(Grads_def(1,:),'o','DisplayName','Finite Diff, default h'); 
plot(Grads_man(1,:),'x','DisplayName','Finite Diff, manual h'); 
plot(df1(eval_points),'-','DisplayName','true val')
legend show
grid on

figure
hold on; 
plot(Grads_def(2,:),'-o','DisplayName','Finite Diff, default h'); 
plot(Grads_man(2,:),'x','DisplayName','Finite Diff, manual h'); 
plot(df2(eval_points),'-','DisplayName','true val')
legend show
grid on

% error
clc
max(abs((Grads_def(1,:) - df1(eval_points))./df1(eval_points)))
max(abs((Grads_man(1,:) - df1(eval_points))./df1(eval_points)))


%% a function to compute Hessians of a function (by finite differences) is also available, see hessian_sparse_grid


%% PART 7: SAVE SPARSE GRID ON FILE

clc

N=3;

aa=[4 1 -2];
bb=[6 5 -1];
knots1=@(n) knots_CC(n,aa(1),bb(1),'nonprob');
knots2=@(n) knots_CC(n,aa(2),bb(2),'nonprob');
knots3=@(n) knots_uniform(n,aa(3),bb(3),'nonprob');
w = 2;
S = create_sparse_grid(N,w,{knots1,knots2,knots3},@lev2knots_doubling);
Sr = reduce_sparse_grid(S);

% save points to 'points.dat'. The first row actually contains two integer
% values, i.e., Sr.size and N
export_sparse_grid_to_file(Sr);

% save points to 'mygrid.dat'
export_sparse_grid_to_file(Sr,'mygrid.dat');

% save points and to 'mygrid_with_weights.dat'
export_sparse_grid_to_file(Sr,'mygrid_with_weights.dat','with_weights');
