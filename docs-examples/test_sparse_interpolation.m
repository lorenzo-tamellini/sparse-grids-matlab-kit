clear

% dimension 
N=6;

% function to be interpolated. input column points, out row vector

f=@(x) 1./(1+0.5/N*sum(x));
% f=@(x) sum(abs(x).^3);


% define sparse grid
[lev2knots,idxset]=define_functions_for_rule('TD',N);
knots=@(n) knots_uniform(n,-1,1,'nonprob');
a=-1; b=1;



% for loop

w_max=6;
interp_error=[];
work=[];


non_grid_points=rand(2000,N)*(b-a)+a;


for w=0:w_max

    disp(w)

    % create grid
    [S,C]=smolyak_grid(N,w,knots,lev2knots,idxset);
    
    Sr=reduce_sparse_grid(S);

    
    % move points to actual interval here point are columns
    % Sr.knots=interval_map(Sr.knots);

    % compute work
    work(w+1)=size(Sr.knots,2);
    
    % compute estimate of polynomial size
    pol_size(w+1)=size(C,1);
    
    % compute the nodal values to be used to interpolate
    function_on_grid=f(Sr.knots);   
    
    % compute interpolated values. Here f_values is column
    f_values = interpolate_on_sparse_grid(S,[],Sr,function_on_grid,non_grid_points);

    % compute error
    interp_error(w+1)=max( abs(  ( f(non_grid_points')' - f_values ) )  ) ;    

end


figure
semilogy(0:w,interp_error,'-or','DisplayName','Numerical Error, w.r.t. grid level');
hold on
semilogy(0:w_max,1./exp(0:w_max),'--o','DisplayName','exp(-level)')
legend show

%%
figure
loglog(work,interp_error,'-or','DisplayName','Numerical Error, w.r.t. #points');
legend show

