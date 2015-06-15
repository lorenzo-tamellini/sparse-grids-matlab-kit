%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------

%% a more complex example of recycle

clc
clear

fs=@(x) sum(x);

N=2; 


previous_evals=[];
Sr_old=[];

for w=0:4
    
    S=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
    Sr= reduce_sparse_grid(S);

    evals_rec = evaluate_on_sparse_grid(fs,Sr,previous_evals,Sr_old);
    
    evals_non_rec = evaluate_on_sparse_grid(fs,Sr);
    
    previous_evals=evals_rec;
    Sr_old=Sr;

    max(abs(evals_non_rec(:)-evals_rec(:))) 

end



%% it is also possible to fine-tune the tolerance for 2 points two be considered equal and hence recyclable

clc
clear

f=@(x) sum(x);

N=2; w=3;
S=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=4;
T=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);


tol=1e-16;
evals_old=evaluate_on_sparse_grid(f,Sr);
evals_r_par =evaluate_on_sparse_grid(f,Tr,evals_old,Sr,1,tol);


%% an extreme test: we recycle from a larger grid to a smaller

clc
clear

f=@(x) sum(x);

N=2; w=3;
S=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Sr= reduce_sparse_grid(S);

w=4;
T=smolyak_grid(N,w,@(n) knots_uniform(n,-1,1),@lev2knots_lin);
Tr= reduce_sparse_grid(T);

evals_nr=evaluate_on_sparse_grid(f,Sr);
evals_r=evaluate_on_sparse_grid(f,Sr,evaluate_on_sparse_grid(f,Tr),Tr);

%[i,j]=max(abs(evals_nr(:)-evals_r(:))) 
max(abs(evals_nr(:)-evals_r(:))) 