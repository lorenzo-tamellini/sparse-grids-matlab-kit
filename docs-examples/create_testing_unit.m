
% .........................................................................
% %% 1 - salvare diversi file per i diversi blocchi di funzioni che testiamo
% insieme? Ora si chiamano tutti 'test_unit.dat'
% %% 2 - Folder tools - NO TEST per converter_functions,
% parallel_toolbox_interface_functions, tye_and_property_check_functions
% .........................................................................

%% test lev2knots_functions (tools/lev2knots_functions)
 
clc
clear
testing_mode = true;

ii = [1 2 3 4];
m_lin = lev2knots_lin(ii);
m_2step = lev2knots_2step(ii);
m_doub = lev2knots_doubling(ii);
m_trip = lev2knots_tripling(ii);
m_GK = lev2knots_GK(ii);


if ~testing_mode
    save('test_unit','m_lin','m_2step','m_doub','m_trip','m_GK')
else
    disp('testing lev2knots')
    L = struct('m_lin',m_lin,'m_2step',m_2step,'m_doub',m_doub,'m_trip',m_trip,'m_GK',m_GK);
    S = load('test_unit','m_lin','m_2step','m_doub','m_trip','m_GK');
    % I could just use isequal(L,S) but this would not give info on which field contains the error, so I have to
    % write up a function for this
    if isequal_sgmk(L,S)
        disp('test on lev2knots function passed')
    end    
end


%% test point generation (tools/knots_functions)
 
% .........................................................................
% %% 2 - cosa facciamo con i triangular leja? Al momento non ci sono ne' nel
% tutorial ne' nella testing uni
% %% 3 - controllare tutorial "equispaced (trapezoidal quadrature rule) and
% midpoint", simboli nel plot fatti di proposito uguali?
% %% 4 - ho aggiunto leja per exp, gamma e beta nel tutorial
% .........................................................................

clc
clear
testing_mode = true;

n=5; a=1; b=4;
[x_unif,w_unif]         = knots_uniform(n,a,b);
[x_CC,w_CC]             = knots_CC(n,a,b);
[x_leja,w_leja]         = knots_leja(n,a,b,'line');
[x_sym_leja,w_sym_leja] = knots_leja(n,a,b,'sym_line');
[x_p_leja,w_p_leja]     = knots_leja(n,a,b,'p_disk');
[x_midp,w_midp]         = knots_midpoint(n,a,b);
[x_trap,w_trap]         = knots_trap(n,a,b);


n=9; mu=0; sig=1;
[x_norm,w_norm]                     = knots_normal(n,mu,sig);
[x_GK,w_GK]                         = knots_GK(n,mu,sig);
[x_norm_Leja,w_norm_Leja]           = knots_normal_leja(n,mu,sig,'line');
[x_norm_sym_Leja,w_norm_sym_Leja]   = knots_normal_leja(n,mu,sig,'sym_line');


n=12; lambda=1; 
[x_exp, w_exp]           = knots_exponential(n,lambda);
[x_exp_leja, w_exp_leja] = knots_exponential(n,lambda);


n=12; alpha=1; beta=2;
[x_gamma, w_gamma]           = knots_gamma(n,alpha,beta);
[x_gamma_leja, w_gamma_leja] = knots_gamma(n,alpha,beta);


n=12; x_a=1; x_b=3; alpha=-0.5; beta=0.5; 
[x_beta,w_beta]                   =knots_beta(n,alpha,beta,x_a,x_b);
[x_beta_leja,w_beta_leja]         = knots_beta_leja(n,alpha,beta,x_a,x_b,'line');
[x_beta_sym_leja,w_beta_sym_leja] = knots_beta_leja(n,alpha,beta,x_a,x_b,'sym_line');


if ~testing_mode
    save('test_unit',...
        'x_unif','w_unif','x_CC','w_CC','x_leja','w_leja','x_sym_leja','w_sym_leja','x_p_leja','w_p_leja','x_midp','w_midp','x_trap','w_trap',...
        'x_norm','w_norm','x_GK','w_GK','x_norm_Leja','w_norm_Leja','x_norm_sym_Leja','w_norm_sym_Leja',...
        'x_exp','w_exp','x_exp_leja','w_exp_leja',...
        'x_gamma','w_gamma','x_gamma_leja','w_gamma_leja',...
        'x_beta','w_beta','x_beta_leja','w_beta_leja','x_beta_sym_leja','w_beta_sym_leja','-append');
else
    disp('testing knots')
    L = struct('x_unif',x_unif,'w_unif',w_unif,'x_CC',x_CC,'w_CC',w_CC,'x_leja',x_leja,'w_leja',w_leja,...
               'x_sym_leja',x_sym_leja,'w_sym_leja',w_sym_leja,'x_p_leja',x_p_leja,'w_p_leja',w_p_leja,...
               'x_midp',x_midp,'w_midp',w_midp,'x_trap',x_trap,'w_trap',w_trap,...
               'x_norm',x_norm,'w_norm',w_norm,'x_GK',x_GK,'w_GK',w_GK,'x_norm_Leja',x_norm_Leja,'w_norm_Leja',w_norm_Leja,'x_norm_sym_Leja',x_norm_sym_Leja,'w_norm_sym_Leja',w_norm_sym_Leja, ...
               'x_exp',x_exp,'w_exp',w_exp,'x_exp_leja',x_exp_leja,'w_exp_leja',w_exp_leja,...
               'x_gamma',x_gamma,'w_gamma',w_gamma,'x_gamma_leja',x_gamma_leja,'w_gamma_leja',w_gamma_leja,...
               'x_beta',x_beta,'w_beta',w_beta,'x_beta_leja',x_beta_leja,'w_beta_leja',w_beta_leja,'x_beta_sym_leja',x_beta_sym_leja,'w_beta_sym_leja',w_beta_sym_leja);
    S = load('test_unit',...
              'x_unif','w_unif','x_CC','w_CC','x_leja','w_leja','x_sym_leja','w_sym_leja','x_p_leja','w_p_leja','x_midp','w_midp','x_trap','w_trap',...
              'x_norm','w_norm','x_GK','w_GK','x_norm_Leja','w_norm_Leja','x_norm_sym_Leja','w_norm_sym_Leja',...
              'x_exp','w_exp','x_exp_leja','w_exp_leja',...
              'x_gamma','w_gamma','x_gamma_leja','w_gamma_leja',...
              'x_beta','w_beta','x_beta_leja','w_beta_leja','x_beta_sym_leja','w_beta_sym_leja');   
    if isequal_sgmk(L,S)
        disp('test on knots passed')
    end    
end


%% test multi-index set generation (tools/idxset_functions)

% .........................................................................
% %% 5 - folder idxset_functions: NO TEST per check_index_admissibility,
% check_set_admissibility, plot_multiidx_set
% %% 6 - non testo define_functions_for_rule ma la uso per generare gli
% input per multiidx_gen. Va bene?
% %% 7 - [lev2nodes,idxset] = define_functions_for_rule(rule,input2). Non
% sarebbe meglio chiamare il secondo output "rule"?
% .........................................................................

clc
clear
testing_mode = true;

jj = [2 3]; min_idx = 0; 
multi_idx_box = multiidx_box_set(jj,min_idx);

N = 2;  
[~,rule_TP]       = define_functions_for_rule('TP',N);
[~,rule_TD]       = define_functions_for_rule('TD',N);
[~,rule_HC]       = define_functions_for_rule('HC',N);
[~,rule_SM]       = define_functions_for_rule('SM',N);

base = 1; w = 3; 
multi_idx_TP      = multiidx_gen(N,rule_TP,w,base); 
multi_idx_TD      = multiidx_gen(N,rule_TD,w,base); 
multi_idx_HC      = multiidx_gen(N,rule_HC,w,base); 
multi_idx_SM      = multiidx_gen(N,rule_SM,w,base);

multi_idx_TD_fast = fast_TD_set(N,w); 

if ~testing_mode
    save('test_unit','multi_idx_box','multi_idx_TP','multi_idx_TD',...
         'multi_idx_HC','multi_idx_SM','multi_idx_TD_fast','-append');
else
    disp('testing multi-index sets')
    L = struct('multi_idx_box',multi_idx_box,'multi_idx_TP',multi_idx_TP,....
               'multi_idx_TD',multi_idx_TD, 'multi_idx_HC',multi_idx_HC,...
               'multi_idx_SM',multi_idx_SM, 'multi_idx_TD_fast',multi_idx_TD_fast);
    S = load('test_unit',...
             'multi_idx_box','multi_idx_TP','multi_idx_TD','multi_idx_HC','multi_idx_SM','multi_idx_TD_fast');   
    if isequal_sgmk(L,S)
        disp('test on multi-index sets passed')
    end    
end


%% test polynomials (tools/polynomials_functions)

% .........................................................................
% %% 1 - vogliamo confrontare le valutazioni in un set di punti, giusto?
% %% 2 - non testo le funzioni per i polinomi di lagrange, ma testo univariate_interpolant
% %% 3 - non abbiamo i polinomi nel tutorial, ma non mi sembra grave
% perche' tanto vengono chiamati solo internamente dalla conversion e sobol
% %% 4 - ok la scelta dei punti in cui valutare i polinomi per lagu and
% hermite? lagu: ln(0.05/(-lambda) ma gen_lagu non so come farlo, tiro a
% occhio?
% %% 5 - quanti punti di valutazione vogliamo mettere?
% .........................................................................

clc
clear
testing_mode = true;

% Legendre and Chebyshev
a = -2; b = 1; n = 20; 
x = linspace(a,b,n); k = 5;  
lege_vals = lege_eval(x,k,a,b);
cheb_vals = cheb_eval(x,k,a,b);

aa = [-2,1]; bb = [1,2]; nn = [20,10];  
[x1,x2] = meshgrid(linspace(aa(1),bb(1),nn(1)),linspace(aa(2),bb(2),nn(2))); 
X = [x1(:),x2(:)]'; kk = [3,5];   
multidim_lege_vals = lege_eval_multidim(X,kk,aa,bb); 
multidim_cheb_vals = cheb_eval_multidim(X,kk,aa,bb);

% Hermite
mi = 0; sigma = 1; n = 50; 
x = linspace(mi-3*sigma,mi+3*sigma,n); k = 5;  
herm_vals = herm_eval(x,k,mi,sigma); 

mmi = [0,1]; ssigma = [1,0.5]; nn = [50,50];  
[x1,x2] = meshgrid(linspace(mmi(1)-3*ssigma(1),mmi(1)+3*ssigma(1),nn(1)),linspace(mmi(2)-3*ssigma(2),mmi(2)+3*ssigma(2),nn(2)));
X = [x1(:),x2(:)]'; kk = [3,5];  
multidim_herm_vals = herm_eval_multidim(X,kk,mmi,ssigma); 

% Laguerre
lambda = 1; n = 20; 
x = linspace(0,log(0.05)/(-lambda),n); k = 5; 
lagu_vals = lagu_eval(x,k,lambda); 

llambda = [1,2]; nn = [20,40];  
[x1,x2] = meshgrid(linspace(0,log(0.05)/(-llambda(1)),nn(1)),linspace(0,log(0.05)/(-llambda(2)),nn(2))); 
X = [x1(:),x2(:)]'; kk = [3,5];  
multidim_lagu_vals = lagu_eval_multidim(X,kk,llambda);

% Generalized Laguerre
alpha = 1; beta = 2; n = 20; 
x = linspace(0,5,n); k = 5; 
generalized_lagu_vals = generalized_lagu_eval(x,k,alpha,beta); 

aalpha = [1,2]; bbeta = [1,1]; nn = [20,40];  
[x1,x2] = meshgrid(linspace(0,5,nn(1)),linspace(0,10,nn(2))); 
X = [x1(:),x2(:)]'; kk = [3,5];  
multidim_generalized_lagu_vals = generalized_lagu_eval_multidim(X,kk,aalpha,bbeta);

% Jacobi
a = -2; b = 1; n = 20; alpha = 0.5; beta = 0.5; 
x = linspace(a,b,n); k = 5;  
jac_vals = jacobi_prob_eval(x,k,alpha,beta,a,b); 

aa = [-2,1]; bb = [1,2]; nn = [20,10]; aalpha = [-0.5,1]; bbeta = [-0.5,2];   
[x1,x2] = meshgrid(linspace(aa(1),bb(1),nn(1)),linspace(aa(2),bb(2),nn(2))); 
X = [x1(:),x2(:)]'; kk = [3,5];   
multidim_jac_vals = jacobi_prob_eval_multidim(X,kk,aalpha,bbeta,aa,bb); 

% univariate interpolant - Lagrange basis
a = 0; b = 2*pi; n = 20; 
x_interp = linspace(a,b,n); % interpolation points
f = @(x) sin(x); 
f_interp = f(x_interp); 
x_eval = linspace(a,b,2*n); % to not take interpolation points only
lagr_vals = univariate_interpolant(x_interp,f_interp,x_eval);  

if ~testing_mode
    save('test_unit',...
         'lege_vals','cheb_vals','multidim_lege_vals','multidim_cheb_vals',...
         'herm_vals','multidim_herm_vals',...
         'lagu_vals','multidim_lagu_vals',...
         'generalized_lagu_vals','multidim_generalized_lagu_vals',...
         'jac_vals','multidim_jac_vals',...
         'lagr_vals','-append');
else
    disp('testing polynomials')
    L = struct('lege_vals',lege_vals,'cheb_vals',cheb_vals,'multidim_lege_vals',multidim_lege_vals,'multidim_cheb_vals',multidim_cheb_vals,...
               'herm_vals',herm_vals,'multidim_herm_vals',multidim_herm_vals,...
               'lagu_vals',lagu_vals,'multidim_lagu_vals',multidim_lagu_vals,...
               'generalized_lagu_vals',generalized_lagu_vals,'multidim_generalized_lagu_vals',multidim_generalized_lagu_vals,...
               'jac_vals',jac_vals,'multidim_jac_vals',multidim_jac_vals,...
               'lagr_vals',lagr_vals);
    S = load('test_unit',...
             'lege_vals','cheb_vals','multidim_lege_vals','multidim_cheb_vals',...
             'herm_vals','multidim_herm_vals',...
             'lagu_vals','multidim_lagu_vals',...
             'generalized_lagu_vals','multidim_generalized_lagu_vals',...
             'jac_vals','multidim_jac_vals', ...
             'lagr_vals');   
    if isequal_sgmk(L,S)
        disp('test on polynomials passed')
    end    
end


%% test sparse grid generation and reduction (main)

% .........................................................................
% %% 1 - NO TEST per plot_sparse_grid,plot_sparse_grid_interpolant,plot3_sparse_grid
% %% 2 - perche' coeff_smolyak = combination_technique(I_smolyak) non si
% puo' inglobare in smolyak_grid_add_multiidx?
% .........................................................................

clc
clear
testing_mode = true;

% given the multi-index set 
I = [
    1 1;
    1 2;
    2 1;
    3 1;
];
knots1 = @(n) knots_uniform(n,0,1);
knots2 = @(n) knots_leja(n,-1,1,'line');
lev2knots = @lev2knots_lin; 
S_given_multiidx = smolyak_grid_multiidx_set(I,{knots1,knots2},lev2knots); 

% given the rule 
N=2; w=3;
knots = @(n) knots_CC(n,-1,1); 
[lev2knots,rule] = define_functions_for_rule('SM',N); 
[S_smolyak,I_smolyak] = smolyak_grid(N,w,knots,lev2knots,rule);
Sr_smolyak = reduce_sparse_grid(S_smolyak); 

% adding one multi-index to S_smolyak
new_idx = [5 1];
coeff_smolyak = combination_technique(I_smolyak); 
[S_add,I_add,coeff_add] = smolyak_grid_add_multiidx(new_idx,S_smolyak,I_smolyak,coeff_smolyak,knots,lev2knots); 

% quick preset
N = 2; w = 3;
[S_quick,Sr_quick] = smolyak_grid_quick_preset(N,w);

if ~testing_mode
    save('test_unit',...
         'S_given_multiidx',...
         'S_smolyak','I_smolyak','Sr_smolyak',...
         'S_add','I_add','coeff_add',...
         'S_quick','Sr_quick','-append');
else
    disp('testing sparse grid generation')
    L = struct('S_given_multiidx',S_given_multiidx,...
               'S_smolyak',S_smolyak,'I_smolyak',I_smolyak,'Sr_smolyak',Sr_smolyak,...
               'S_add',S_add,'I_add',I_add,'coeff_add',coeff_add,...
               'S_quick',S_quick,'Sr_quick',Sr_quick);
    S = load('test_unit',...
             'S_given_multiidx',...
             'S_smolyak','I_smolyak','Sr_smolyak',...
             'S_add','I_add','coeff_add', ...
             'S_quick','Sr_quick');   
    if isequal_sgmk(L,S)
        disp('test on sparse grid generation passed')
    end    
end

%% test on quadrature

clear





%%%%%%%%%%%% USE KEEP.M AND GIVE CREDIT  %%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iseq = isequal_sgmk(L,S)

    tested_fields = fieldnames(S);
    nb_tests = length(tested_fields);
    for t = 1:nb_tests
        if ~isequal(L.(tested_fields{t}), S.(tested_fields{t}) )
            iseq = false;
            error(['error in test of ',tested_fields{t}])
        end
    end
    iseq = true;
    
end
