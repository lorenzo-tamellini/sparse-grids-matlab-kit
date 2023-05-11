
% RUN_TESTING_UNIT  
% Run the testing unit that compares function output with the corresponding
% results stored in the following files
% 
% test_unit_lev2knots.mat
% test_unit_knots.mat
% test_unit_multiidx_set.mat
% test_unit_polynomials.mat
% test_unit_grid_gen_and_red.mat
% test_unit_evaluate.mat
% test_unit_quadrature.mat
% test_unit_interpolate.mat
% test_unit_gPCE.mat
% test_unit_sobol.mat
% test_unit_gradient_and_hessian.mat
% test_unit_adaptive.mat

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------

global MATLAB_SPARSE_KIT_VERBOSE
MATLAB_SPARSE_KIT_VERBOSE=0;

clc

%% test lev2knots_functions (tools/lev2knots_functions)
 
clear

ii = [1 2 3 4];
m_lin = lev2knots_lin(ii);
m_2step = lev2knots_2step(ii);
m_doub = lev2knots_doubling(ii);
m_trip = lev2knots_tripling(ii);
m_GK = lev2knots_GK(ii);

disp('== testing lev2knots ==')
L = struct('m_lin',m_lin,'m_2step',m_2step,'m_doub',m_doub,'m_trip',m_trip,'m_GK',m_GK);
S = load('test_unit_lev2knots','m_lin','m_2step','m_doub','m_trip','m_GK');
if isequal_sgmk(L,S) % this function is like isequal between struct but 1) it is robust to eps-machine noise and 2) gives some info on which fields are different
    disp('test on lev2knots function passed')
end    

%% test point generation (tools/knots_functions)

clear

% uniform pdf
n = 5; a = 1; b =4 ;
[x_unif,w_unif]         = knots_uniform(n,a,b);
[x_CC,w_CC]             = knots_CC(n,a,b);
[x_leja,w_leja]         = knots_leja(n,a,b,'line');
[x_sym_leja,w_sym_leja] = knots_leja(n,a,b,'sym_line');
[x_p_leja,w_p_leja]     = knots_leja(n,a,b,'p_disk');
[x_midp,w_midp]         = knots_midpoint(n,a,b);
[x_trap,w_trap]         = knots_trap(n,a,b);

% normal pdf
n = 9; mu = 0; sigma = 1;
[x_norm,w_norm]                     = knots_normal(n,mu,sigma);
[x_GK,w_GK]                         = knots_GK(n,mu,sigma);
[x_norm_Leja,w_norm_Leja]           = knots_normal_leja(n,mu,sigma,'line');
[x_norm_sym_Leja,w_norm_sym_Leja]   = knots_normal_leja(n,mu,sigma,'sym_line');

% exponential pdf
n = 12; lambda = 1; 
[x_exp, w_exp]           = knots_exponential(n,lambda);
[x_exp_leja, w_exp_leja] = knots_exponential(n,lambda);

% gamma pdf
n = 12; alpha = 1; beta = 2;
[x_gamma, w_gamma]           = knots_gamma(n,alpha,beta);
[x_gamma_leja, w_gamma_leja] = knots_gamma(n,alpha,beta);

% beta pdf
n = 12; a = 1; b = 3; 
alpha = -0.5; beta = 0.5; 
[x_beta,w_beta]                   = knots_beta(n,alpha,beta,a,b);
[x_beta_leja,w_beta_leja]         = knots_beta_leja(n,alpha,beta,a,b,'line');
alpha = 1.5; beta = 1.5; % alpha = beta for symmetric points 
[x_beta_sym_leja,w_beta_sym_leja] = knots_beta_leja(n,alpha,beta,a,b,'sym_line');

% triangular pdf
n = 12; a = 0; b = 2;
[x_triang,w_triang] = knots_triangular_leja(n,a,b);

disp('== testing knots ==')
L = struct('x_unif',x_unif,'w_unif',w_unif,'x_CC',x_CC,'w_CC',w_CC,'x_leja',x_leja,'w_leja',w_leja,...
           'x_sym_leja',x_sym_leja,'w_sym_leja',w_sym_leja,'x_p_leja',x_p_leja,'w_p_leja',w_p_leja,...
           'x_midp',x_midp,'w_midp',w_midp,'x_trap',x_trap,'w_trap',w_trap,...
           'x_norm',x_norm,'w_norm',w_norm,'x_GK',x_GK,'w_GK',w_GK,'x_norm_Leja',x_norm_Leja,'w_norm_Leja',w_norm_Leja,...
           'x_norm_sym_Leja',x_norm_sym_Leja,'w_norm_sym_Leja',w_norm_sym_Leja, ...
           'x_exp',x_exp,'w_exp',w_exp,'x_exp_leja',x_exp_leja,'w_exp_leja',w_exp_leja,...
           'x_gamma',x_gamma,'w_gamma',w_gamma,'x_gamma_leja',x_gamma_leja,'w_gamma_leja',w_gamma_leja,...
           'x_beta',x_beta,'w_beta',w_beta,'x_beta_leja',x_beta_leja,'w_beta_leja',w_beta_leja,...
           'x_beta_sym_leja',x_beta_sym_leja,'w_beta_sym_leja',w_beta_sym_leja,...
           'x_triang',x_triang,'w_triang',w_triang);
S = load('test_unit_knots',...
          'x_unif','w_unif','x_CC','w_CC','x_leja','w_leja','x_sym_leja','w_sym_leja','x_p_leja','w_p_leja','x_midp','w_midp','x_trap','w_trap',...
          'x_norm','w_norm','x_GK','w_GK','x_norm_Leja','w_norm_Leja','x_norm_sym_Leja','w_norm_sym_Leja',...
          'x_exp','w_exp','x_exp_leja','w_exp_leja',...
          'x_gamma','w_gamma','x_gamma_leja','w_gamma_leja',...
          'x_beta','w_beta','x_beta_leja','w_beta_leja','x_beta_sym_leja','w_beta_sym_leja',...
          'x_triang','w_triang');   
if isequal_sgmk(L,S)
    disp('test on knots passed')
end    


%% test multi-index set generation (tools/idxset_functions)

clear

jj = [2 3]; min_idx = 0; 
multi_idx_box = multiidx_box_set(jj,min_idx);

N = 2;  
% isotropic sets
[~,rule_TP]       = define_functions_for_rule('TP',N);
[~,rule_TD]       = define_functions_for_rule('TD',N);
[~,rule_HC]       = define_functions_for_rule('HC',N);
[~,rule_SM]       = define_functions_for_rule('SM',N);
% anisotropic sets
rates = [2,3]; 
[~,rule_SM_aniso] = define_functions_for_rule('SM',rates);

base = 1; w = 3; 
multi_idx_TP       = multiidx_gen(N,rule_TP,w,base); 
multi_idx_TD       = multiidx_gen(N,rule_TD,w,base); 
multi_idx_HC       = multiidx_gen(N,rule_HC,w,base); 
multi_idx_SM       = multiidx_gen(N,rule_SM,w,base);
multi_idx_SM_aniso = multiidx_gen(N,rule_SM,w,base);

% fast TD
multi_idx_TD_fast = fast_TD_set(N,w); 

disp('== testing multi-index sets ==')
L = struct('multi_idx_box',multi_idx_box,'multi_idx_TP',multi_idx_TP,....
           'multi_idx_TD',multi_idx_TD,'multi_idx_HC',multi_idx_HC,...
           'multi_idx_SM',multi_idx_SM,'multi_idx_SM_aniso',multi_idx_SM_aniso,...
           'multi_idx_TD_fast',multi_idx_TD_fast);
S = load('test_unit_multiidx_set','multi_idx_box','multi_idx_TP','multi_idx_TD',...
         'multi_idx_HC','multi_idx_SM','multi_idx_SM_aniso','multi_idx_TD_fast');   
if isequal_sgmk(L,S)
    disp('test on multi-index sets passed')
end    

%% test polynomials (tools/polynomials_functions)

clear

k = 5; kk = [3,5]; % polynomial order (univariate and bi-variate case, respectively)
load('test_unit_polynomials','x_unif','X_unif_multidim','x_norm','X_norm_multidim',...
     'x_exp','X_exp_multidim','x_gamma','X_gamma_multidim','x_beta','X_beta_multidim',...
     'x_interp','x_eval'); 

% Legendre and Chebyshev
a = -2; b = 1; 
lege_vals = lege_eval(x_unif,k,a,b);
cheb_vals = cheb_eval(x_unif,k,a,b);

aa = [-2,1]; bb = [1,2]; 
multidim_lege_vals = lege_eval_multidim(X_unif_multidim,kk,aa,bb); 
multidim_cheb_vals = cheb_eval_multidim(X_unif_multidim,kk,aa,bb);

% Hermite
mi = 0; sigma = 1;  
herm_vals = herm_eval(x_norm,k,mi,sigma); 

mmi = [0,1]; ssigma = [1,0.5]; 
multidim_herm_vals = herm_eval_multidim(X_norm_multidim,kk,mmi,ssigma); 

% Laguerre
lambda = 2; 
lagu_vals = lagu_eval(x_exp,k,lambda); 

llambda = [1,2]; 
multidim_lagu_vals = lagu_eval_multidim(X_exp_multidim,kk,llambda);

% Generalized Laguerre
alpha = 1; beta = 2; 
generalized_lagu_vals = generalized_lagu_eval(x_gamma,k,alpha,beta); % parameters alpha-1,1/beta for consistency with Matlab definitions

aalpha = [1,2]; bbeta = [1,0.5];  
multidim_generalized_lagu_vals = generalized_lagu_eval_multidim(X_gamma_multidim,kk,aalpha,bbeta);

% Jacobi
a = -2; b = 1; 
alpha = -0.5; beta = 0.5; 
jac_vals = jacobi_prob_eval(x_beta,k,alpha,beta,a,b); 

aa = [-2,1]; bb = [1,2]; 
aalpha = [-0.5,1]; bbeta = [-0.5,2];   
multidim_jac_vals = jacobi_prob_eval_multidim(X_beta_multidim,kk,aalpha,bbeta,aa,bb); 

% univariate interpolant - Lagrange basis
f = @(x) sin(x); 
f_interp = f(x_interp); 
lagr_vals = univariate_interpolant(x_interp,f_interp,x_eval);  

disp('== testing polynomials ==')
L = struct('lege_vals',lege_vals,'cheb_vals',cheb_vals,'multidim_lege_vals',multidim_lege_vals,'multidim_cheb_vals',multidim_cheb_vals,...
           'herm_vals',herm_vals,'multidim_herm_vals',multidim_herm_vals,...
           'lagu_vals',lagu_vals,'multidim_lagu_vals',multidim_lagu_vals,...
           'generalized_lagu_vals',generalized_lagu_vals,'multidim_generalized_lagu_vals',multidim_generalized_lagu_vals,...
           'jac_vals',jac_vals,'multidim_jac_vals',multidim_jac_vals,...
           'lagr_vals',lagr_vals);
S = load('test_unit_polynomials',...
         'lege_vals','cheb_vals','multidim_lege_vals','multidim_cheb_vals',...
         'herm_vals','multidim_herm_vals',...
         'lagu_vals','multidim_lagu_vals',...
         'generalized_lagu_vals','multidim_generalized_lagu_vals',...
         'jac_vals','multidim_jac_vals', ...
         'lagr_vals');   
if isequal_sgmk(L,S)
    disp('test on polynomials passed')
end    

%% test sparse grid generation and reduction (main)

clear

% given the multi-index set 
I = [
    1 1;
    1 2;
    2 1;
    3 1;
];
knots1           = @(n) knots_uniform(n,0,1);
knots2           = @(n) knots_leja(n,-1,1,'line');
lev2knots        = @lev2knots_lin; 
S_given_multiidx = create_sparse_grid_multiidx_set(I,{knots1,knots2},lev2knots); 

% given the rule 
N=2; w=3;
knots                 = @(n) knots_CC(n,-1,1); 
[lev2knots,rule]      = define_functions_for_rule('SM',N); 
[S_smolyak,I_smolyak] = create_sparse_grid(N,w,knots,lev2knots,rule);
Sr_smolyak            = reduce_sparse_grid(S_smolyak); 

% adding one multi-index to S_smolyak
new_idx = [5 1];
coeff_smolyak           = combination_technique(I_smolyak); 
[S_add,I_add,coeff_add] = create_sparse_grid_add_multiidx(new_idx,S_smolyak,I_smolyak,coeff_smolyak,knots,lev2knots); 

% quick preset
N = 2; w = 3;
[S_quick,Sr_quick] = create_sparse_grid_quick_preset(N,w);

disp('== testing sparse grid generation and reduction ==')
L = struct('S_given_multiidx',S_given_multiidx,...
           'S_smolyak',S_smolyak,'I_smolyak',I_smolyak,'Sr_smolyak',Sr_smolyak,...
           'S_add',S_add,'I_add',I_add,'coeff_add',coeff_add,...
           'S_quick',S_quick,'Sr_quick',Sr_quick);
S = load('test_unit_grid_gen_and_red',...
         'S_given_multiidx',...
         'S_smolyak','I_smolyak','Sr_smolyak',...
         'S_add','I_add','coeff_add', ...
         'S_quick','Sr_quick');   
     

if isequal_sgmk(L,S)
    disp('test on sparse grid generation and reduction passed')
end    

%% test on function evaluation on sparse grid (main)

clear

f = @(x) sum(x);

N = 2; w = 3;
S  = create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr = reduce_sparse_grid(S);

f_evals = evaluate_on_sparse_grid(f,Sr);

w = 4;
T  = create_sparse_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr = reduce_sparse_grid(T);

f_evals_rec = evaluate_on_sparse_grid(f,T,Tr,f_evals,S,Sr);

disp('== testing function evaluation on sparse grid ==')
L = struct('f_evals',f_evals,'f_evals_rec',f_evals_rec);
S = load('test_unit_evaluate','f_evals','f_evals_rec');   
if isequal_sgmk(L,S)
    disp('test on function evaluation on sparse grid passed')
end    

%% test on integration on sparse grid (main)

clear

f = @(x) prod(1./sqrt(x+3));

N = 4; w =4 ;
knots = @(n) knots_CC(n,-1,1);
S     = create_sparse_grid(N,w,knots,@lev2knots_doubling);
Sr    = reduce_sparse_grid(S);

f_quad = quadrature_on_sparse_grid(f,Sr); 

disp('== testing integration on sparse grid ==')
L = struct('f_quad',f_quad);
S = load('test_unit_quadrature','f_quad');   
if isequal_sgmk(L,S)
    disp('test on integration on sparse grid passed')
end

%% test on interpolation on sparse grid (main)

clear

f = @(x) prod(1./sqrt(x+3)); 

N = 2; w = 4;
knots = @(n) knots_CC(n,-1,1);
S     = create_sparse_grid(N,w,knots,@lev2knots_doubling);
Sr    = reduce_sparse_grid(S);

x               = linspace(-1,1,10); 
[X1,X2]         = meshgrid(x,x); 
non_grid_points = [X1(:),X2(:)]'; 
f_on_grid       = evaluate_on_sparse_grid(f, Sr);
f_values        = interpolate_on_sparse_grid(S,Sr,f_on_grid,non_grid_points);

disp('== testing interpolation on sparse grid ==')
L = struct('f_values',f_values);
S = load('test_unit_interpolate','f_values');   
if isequal_sgmk(L,S)
    disp('test on interpolation on sparse grid passed')
end    

%% test on generalized Polynomial Chaos Expansion (gPCE) 

clear

f = @(x) prod(1./sqrt(x+3));

N = 2; w = 5; a = -1; b = 1; 
knots     = @(n) knots_uniform(n,a,b); 
lev2knots = @lev2knots_lin; 
idxset    = @(i) prod(i); 
S         = create_sparse_grid(N,w,knots,lev2knots,idxset);
Sr        = reduce_sparse_grid(S);
f_on_grid = evaluate_on_sparse_grid(f,Sr);

[PCE_coeffs,PCE_multiidx] = convert_to_modal(S,Sr,f_on_grid,[a a; b b],'legendre');

disp('== testing gPCE ==')
L = struct('PCE_coeffs',PCE_coeffs,'PCE_multiidx',PCE_multiidx);
S = load('test_unit_gPCE','PCE_coeffs','PCE_multiidx');   
if isequal_sgmk(L,S)
    disp('test on gPCE passed')
end    

%% test on Sobol indices 

clear

f = @(x) 1./(1 + 5*x(1,:).^2 + x(2,:).^2 + x(3,:).^2); 

N = 3; w = 5; a = -1; b = 1; 
knots     = @(n) knots_uniform(n,a,b); 
lev2knots = @lev2knots_lin; 
idxset    = @(i) prod(i); 
S         = create_sparse_grid(N,w,knots,lev2knots,idxset);
Sr        = reduce_sparse_grid(S);
f_on_grid = evaluate_on_sparse_grid(f,Sr);

[Sob_i,Tot_Sob_i,Mean,Var] = compute_sobol_indices_from_sparse_grid(S,Sr,f_on_grid,[a a a; b b b],'legendre');

disp('== testing Sobol indices ==')
L = struct('Sob_i',Sob_i,'Tot_Sob_i',Tot_Sob_i,'Mean',Mean,'Var',Var);
S = load('test_unit_sobol','Sob_i','Tot_Sob_i','Mean','Var');   
if isequal_sgmk(L,S)
    disp('test on Sobol indices passed')
end    

%% test on gradient and Hessian

clear

f = @(x) 1./(1+0.5*sum(x.^2)); 

N = 2; aa = [4 1]; bb = [6 5]; domain = [aa; bb]; 
w = 4;
knots1    = @(n) knots_CC(n,aa(1),bb(1));
knots2    = @(n) knots_CC(n,aa(2),bb(2));
S         = create_sparse_grid(N,w,{knots1,knots2},@lev2knots_doubling);
Sr        = reduce_sparse_grid(S);
f_on_grid = evaluate_on_sparse_grid(f,Sr);

x1        = linspace(aa(1),bb(1),10);
x2        = linspace(aa(2),bb(2),30);
[X1,X2]   = meshgrid(x1,x2); 
eval_points = [X1(:),X2(:)]'; 

grad = derive_sparse_grid(S,Sr,f_on_grid,domain,eval_points);
eval_point_hess = eval_points(:,10); % pick one point
Hess = hessian_sparse_grid(S,Sr,f_on_grid,domain,eval_point_hess);

disp('== testing gradients and derivatives ==')
L = struct('grad',grad,'Hess',Hess);
S = load('test_unit_gradient_and_hessian','grad','Hess');   
if isequal_sgmk(L,S)
    disp('test on gradients and derivatives passed')
end    

%% test on adaptive algorithm

clear

f=@(x) 1./(x(1)^2+x(2)^2 + 0.3);

N = 2; 
a = -1; b = 1;
knots     = @(n) knots_CC(n,a,b);
lev2knots = @lev2knots_doubling;

controls.paral     = NaN;
controls.max_pts   = 200;
controls.prof_toll = 1e-10;
prev_adapt         = [];
controls.nested    = true;

% different profit indicators 
controls.prof    = 'Linf'; 
S_Linf = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

controls.prof    = 'Linf/new_points'; 
S_Linf_newpts = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

controls.prof    = 'deltaint'; 
S_deltaint = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

controls.prof    = 'deltaint/new_points'; 
S_deltaint_newpts = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

controls.pdf = @(Y) prod(1/(b-a) * ones(size(Y))); 
controls.prof    = 'weighted Linf'; 
S_weighLinf = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);

controls.prof    = 'weighted Linf/new_points'; 
S_weighLinf_newpts = adapt_sparse_grid(f,N,knots,lev2knots,prev_adapt,controls);


% buffering 
N = 3; 
f = @(x) prod(x);
knots     = @(n) knots_CC(n, 0, 2);
lev2knots = @lev2knots_doubling;
controls.nested = true;

controls.paral     = NaN; 
controls.max_pts   = 200;
controls.prof_toll = 1e-10;
prev_adapt         = [];
controls.plot      = false;

controls.var_buffer_size = 2;
S_buff = adapt_sparse_grid(f,N,knots, lev2knots, prev_adapt, controls);

disp('== testing adaptive algorithm for sparse grid generation ==')
L = struct('S_Linf',S_Linf,'S_Linf_newpts',S_Linf_newpts,...
           'S_weighLinf',S_weighLinf,'S_weighLinf_newpts',S_weighLinf_newpts,...
           'S_deltaint',S_deltaint,'S_deltaint_newpts',S_deltaint_newpts,...
           'S_buff',S_buff);
S = load('test_unit_adaptive','S_Linf','S_Linf_newpts',...
         'S_weighLinf','S_weighLinf_newpts',...
         'S_deltaint','S_deltaint_newpts','S_buff');   
     
if isequal_sgmk(L,S)
    disp('test on adaptive algorithm for sparse grid generation passed')
end    


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iseq = isequal_sgmk(L,S)

    tested_fields = fieldnames(S);
    nb_tests = length(tested_fields);
    
    for t = 1:nb_tests

        if isstruct(L.(tested_fields{t})) % unpack test on each field with recursive call 
            len = length(L.(tested_fields{t})); % the field could be a struct array, not just a struct
            for s = 1:len
                iseq = isequal_sgmk(L.(tested_fields{t})(s),S.(tested_fields{t})(s));
            end
            
        elseif isreal(L.(tested_fields{t})) % test up to tol
            if max( abs( L.(tested_fields{t})-S.(tested_fields{t}) ) ) > 1e-15 
                iseq = false;
                error(['error in test of ',tested_fields{t}])
            end
            
        elseif iscell(L.(tested_fields{t})) % this is e.g. knots per dim
            % convert to number and compare
            Lcell2mat = cell2mat(L.(tested_fields{t}));
            Scell2mat = cell2mat(S.(tested_fields{t}));
            if max(abs( Lcell2mat-Scell2mat ) ) > 1e-15 
                iseq = false;
                error(['error in test of ',tested_fields{t}])
            end
            
        else
            error(['the field ',tested_fields{t},' is neither struct nor number, fix this'])
        end
    end
    iseq = true;
    
end