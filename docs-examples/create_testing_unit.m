

%% test lev2knots_functions
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

%% test point generation


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


if ~testing_mode
    save('test_unit',...
        'x_unif','w_unif','x_CC','w_CC','x_leja','w_leja','x_sym_leja','w_sym_leja','x_p_leja','w_p_leja','x_midp','w_midp','x_trap','w_trap',...
        'x_norm','w_norm','x_GK','w_GK','x_norm_Leja','w_norm_Leja','x_norm_sym_Leja','w_norm_sym_Leja','-append');
else
    disp('testing knots')
    L = struct('x_unif',x_unif,'w_unif',w_unif,'x_CC',x_CC,'w_CC',w_CC,'x_leja',x_leja,'w_leja',w_leja,...
               'x_sym_leja',x_sym_leja,'w_sym_leja',w_sym_leja,'x_p_leja',x_p_leja,'w_p_leja',w_p_leja,...
               'x_midp',x_midp,'w_midp',w_midp,'x_trap',x_trap,'w_trap',w_trap,...
               'x_norm',x_norm,'w_norm',w_norm,'x_GK',x_GK,'w_GK',w_GK,'x_norm_Leja',x_norm_Leja,'w_norm_Leja',w_norm_Leja,'x_norm_sym_Leja',x_norm_sym_Leja,'w_norm_sym_Leja',w_norm_sym_Leja);
    S = load('test_unit',...
        'x_unif','w_unif','x_CC','w_CC','x_leja','w_leja','x_sym_leja','w_sym_leja','x_p_leja','w_p_leja','x_midp','w_midp','x_trap','w_trap',...
        'x_norm','w_norm','x_GK','w_GK','x_norm_Leja','w_norm_Leja','x_norm_sym_Leja','w_norm_sym_Leja');   
    if isequal_sgmk(L,S)
        disp('test on knots passed')
    end    
end


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
