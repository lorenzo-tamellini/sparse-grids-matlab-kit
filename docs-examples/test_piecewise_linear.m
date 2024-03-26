warning('NON GRID POINTS VANNO ORDINATI')
warning('I PUNTI DI GRIGLIA SONO ORDINATI?')

%% --------------------------------------------------------------
% test 1d piecewise lin fast 

% generate sorted points for simplicity

A = 0;
B = 1;
grid_pts = 1;

all_knots = linspace(A,B,grid_pts);

non_grid_points = sort(rand(1,1500));
% non_grid_points = linspace(A,B,100)';

central_knot = 1;
L = piecewise_lin_eval_fast(central_knot,all_knots,non_grid_points,length(non_grid_points));

plot(non_grid_points,L,'-x')
hold on
plot(all_knots,0*all_knots,'o','LineWidth',2)

%% faster to sort column or row vector? they are more or less identical, and very fast anyway

disp('------------------------')
vec = rand(1,100000);
vec_col = vec';

time_row = [];
time_col = [];
for k = 1:10
    tic;
    sort(vec);
    t = toc;
    time_row(end+1) = t;
    tic;
    sort(vec_col);
    t = toc;
    time_col(end+1) = t;
    
end
mean(time_row)
mean(time_col)



%% finally ready

f = @(y) exp(2*y(1) + 2*y(2));

N = 2;
w = 3;
knots = @(n) knots_trap(n,0,1);
lev2knots= @lev2knots_doubling;

S = create_sparse_grid(N,w,knots,lev2knots);
Sr = reduce_sparse_grid(S);
function_on_Sr = evaluate_on_sparse_grid(f,Sr);

% plot_sparse_grid(Sr)

% % first a cloud of random points
% non_grid_points = rand(2,1000);
% f_interp_on_rand_set = interpolate_on_sparse_grid_piecewise_lin(S,Sr,function_on_grid,non_grid_points); 
% plot3(non_grid_points(1,:),non_grid_points(2,:),f_interp_on_rand_set,'xr')
% hold on

% then a struct mesh for surf
NP=20;
xp = linspace(0,1,NP);
yp = linspace(0,1,NP);
[XP,YP] = meshgrid(xp,yp);

nb_pts = length(xp)*length(yp);
PTS = zeros(2,nb_pts);
PTS(1,:)= XP(:)';
PTS(2,:)= YP(:)';

% interpolate on sparse grid
f_interp_on_surf_grid = interpolate_on_sparse_grid_piecewise_lin(S,Sr,function_on_Sr,PTS); 

% reshape to use surf
FIP = reshape(f_interp_on_surf_grid,size(XP));
surf(XP,YP,FIP);
xlabel('y_1')
ylabel('y_2')
hold on
plot_sparse_grid(Sr)
plot3(Sr.knots(1,:),Sr.knots(2,:),function_on_Sr,'ok','MarkerFaceColor','k')



%% convergence compared to a lagrangian grid

clear
f = @(y) exp(y(1) + y(2));
%f = @(y) 0.*(y(1)<=0.5 && y(2)<=0.5) + 1.*(y(1)>0.5 && y(2)<=0.5) + 2.*(y(1)<=0.5 && y(2)>0.5) + 3.*(y(1)>0.5 && y(2)>0.5);
%f = @(y) (atan((y(1)-0.45)*60)+3).*(atan((y(2)-0.45)*60)+3);

N = 2;
wmax = 6;
lev2knots= @lev2knots_doubling;

% a struct mesh for errors
NP=50;
xp = linspace(0,1,NP);
yp = linspace(0,1,NP);
[XP,YP] = meshgrid(xp,yp);

nb_pts = length(xp)*length(yp);
PTS = zeros(2,nb_pts);
PTS(1,:)= XP(:)';
PTS(2,:)= YP(:)';

function_on_PTS = evaluate_on_sparse_grid(f,asreduced(PTS));

err_lin=zeros(1,wmax);
err_lagr=zeros(1,wmax);

for w = 1:wmax
    
    disp(w)
    knots_for_lin = @(n) knots_trap(n,0,1);
    knots_for_lagr = @(n) knots_CC(n,0,1);
    
    S_lin = create_sparse_grid(N,w,knots_for_lin,lev2knots);
    Sr_lin = reduce_sparse_grid(S_lin);
    function_on_grid_lin = evaluate_on_sparse_grid(f,Sr_lin);
    
    f_interp_on_grid_lin = interpolate_on_sparse_grid_piecewise_lin(S_lin,Sr_lin,function_on_grid_lin,PTS);
    
    S_lagr = create_sparse_grid(N,w,knots_for_lagr,lev2knots);
    Sr_lagr = reduce_sparse_grid(S_lagr);
    function_on_grid_lagr = evaluate_on_sparse_grid(f,Sr_lagr);
    
    f_interp_on_grid_lagr = interpolate_on_sparse_grid(S_lagr,Sr_lagr,function_on_grid_lagr,PTS);
    
    err_lin(w) = max(abs(f_interp_on_grid_lin - function_on_PTS)); 
    err_lagr(w) = max(abs(f_interp_on_grid_lagr - function_on_PTS));

end


semilogy(1:wmax,err_lin,'-o','DisplayName','linear');
hold on
semilogy(1:wmax,err_lagr,'-o','DisplayName','lagr');
legend show

figure
plot3(PTS(1,:),PTS(2,:),f_interp_on_grid_lagr,'xr')
hold on
plot3(PTS(1,:),PTS(2,:),f_interp_on_grid_lin,'xb')


%% 0) marcare il commit attuale come release 23.5 (come citare la versione running?)
% 0) branch
% 0**) sistemare bug attuale
% 1**) sort dei punti per dimensione, ricordarsi il sorter
% 2**) un-sort delle valutazioni prima della moltiplicazione sulle dimensioni
% 3) flag "piecewise" nel create sparse grid, nel reduce, nella as_reduced etc
% 4) usare i flag per avere una sola funzione interpolate
% 5) commenti, testing unit, tutorial
% 6) i punti 1d a volte sono decrescenti ... non funziona per quelli
% 7) piecewise multidim