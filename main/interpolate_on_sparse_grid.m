function f_values = interpolate_on_sparse_grid(S,this_was_once_interval_map,Sr,function_on_grid,non_grid_points) %#ok<INUSL>

% F_VALUES = INTERPOLATE_ON_SPARSE_GRID(S,SR,FUNCTION_ON_GRID,NON_GRID_POINTS) interpolates 
%   a vector-valued function F: R^N -> R^V defined on the sparse grid S. S is a sparse grid 
%   and SR is its reduced counterpart.
%   FUNCTION_ON_GRID is a matrix containing the evaluation of F on the points of SR. 
%   Its dimensions are: (number_of_points_in_the_sparse_grid) X V
%   NON_GRID_POINTS is the set of points (not belonging to the sparse grid) where one wants to 
%   evaluate the function. It is a matrix, each row is a different point.
%   F_VALUES is a matrix containing the evaluation of the vector-valued function F 
%   in each of the non_grid_points. Its dimensions are: (number_of_non_grid_point) x V

if nargin == 5
    errmsg=['too many input arguments. This is likely due to the fact that the new version of '...
        'interpolate_on_sparse_grid does not accept any longer INTERVAL_MAP as second input argument. '...
        'The new function call is '...
        'F_VALUES = INTERPOLATE_ON_SPARSE_GRID(S,SR,FUNCTION_ON_GRID,NON_GRID_POINTS). '...
        'To fix this, either regenerate your sparse grid using INTERVAL_MAP as input (see help SMOLYAK_GRID) '...
        'or modify the fields KNOTS of your sparse grid and reduced sparse grid applying INTERVAL_MAP '...
        '(e.g.  S(i).knots = interval_map(S(i).knots) ) '...
        'This message will disappear in future releases of SPARSE_GRID_MATLAB_KIT'];
    error(errmsg)
end


V = size(function_on_grid,2);
nb_pts   = size(non_grid_points,1);
f_values = zeros( nb_pts , V);

nb_grids=length(S);

% I need to go back from each point of each grid to its corresponding value in function_on_grid (F comes from an evaluation over a reduced grid)
% I only have a global mapping from [S.knots] to Sr.knots, so I need a global counter that scrolls [S.knots]
global_knot_counter=1;

% loop over the grids
for i=1:nb_grids
    
    % some of the grids in S structure are empty, so I skip it
    if isempty(S(i).weights)
        continue
    end
    
    % this is the set of points where I build the tensor lagrange function
    knots=S(i).knots;
    if ~isempty(interval_map)
        knots=interval_map( knots );
    end
    
    % I will need the knots in each dimension separately, to collocate the lagrange function.
    % I compute them once for all. As the number of knots is
    % different in each direction, I use a cell array
    dim_tot=size(knots,1);
    knots_per_dim=cell(1,dim_tot);
    
    for dim=1:dim_tot
        knots_per_dim{dim}=unique(knots(dim,:));
    end
    
    
    % I also have to take into account the coefficient of the sparse grid
    
    coeff=S(i).coeff;
    
    
    % we could just loop on the interpolation points of the current tensor grid, 
    % and for each one evaluate the corresponding lagrangian polynomial on
    % the set of non_grid_points, but this is not convenient. Actually doing this we recompute the same thing over and over! 
    % think e.g. of the lagrangian functions for the knot [x1 y1 z1] and
    % for the knot [x1 y1 z2]: the parts on x and y are the same!
    
    % Therefore we evaluate each monodim_lagr just once and save all these
    % evaluations in a cell array. We will combine them suitably afterward.
    %
    % Such cell array stores one matrix per direction, i.e. N matrices. 
    % In turn, the n-th matrix contains the evaluations of each monodimensional lagrange polynomial 
    % in direction n on all the n-th coordinates of the non_grid_points.
    % Its dimension will be therefore (number_of_points_to_evaluate) X (number_of_lagrangian_polynomials_in_the_nth_direction)
    % 
    %
    % This is actually all the information that we will need to combine to
    % get the final interpolant value.
    
    mono_lagr_eval=cell(1,dim_tot);
    
    % loop on directions
    for dim=1:dim_tot
        
        % this is how many grid points in the current direction for the
        % current tensor grid
        K=length(knots_per_dim{dim});
        
        % allocate space for evaluations, conisistently with the shape of
        % non_grid_points (i.e. one point per row, that is all dim-th
        % coordinates in 1 column)
        mono_lagr_eval{dim}=zeros(nb_pts,K);
        
        % loop on each node of the current dimension and evaluate the corresponding monodim lagr polynomial.
        % We will need an auxiliary vector to pick up the current knot (where the lagr pol is centered) and
        % the remaining knots (where the lagr pol is zero)
        aux=1:K;
        for k=aux
            mono_lagr_eval{dim}(:,k) = lagr_eval(knots_per_dim{dim}(k),knots_per_dim{dim}(aux~=k),non_grid_points(:,dim));
        end
        
    end
    
    % now put everything together. We have to take the tensor product of
    % each of the monodim lagr pol we have evaluated. That is, we have to
    % pick one column for each matrix in the cell array and dot-multiply them.
    % all the possible combinations have to be generated !
    %
    % once this is done, we have the evaluation of each multidim lagr
    % polynomial on the non_grid_points, which we will then multiply by the
    % corresponding nodal value and eventually sum everything up.

    % We start by generating the combination of columns. We actually don't
    % need to generate them, but only to recover it from the matrix knots,
    % which already contains all the points of the grid, i.e. all the
    % combinations of 1D points!
    %
    % Given a matrix of points like
    %
    % knots=[a1 a1 b1 b1 a1 a1 b1 b1 ...
    %        a2 b2 a2 b2 .....
    %        
    % combi is
    %
    % combi=[1 1 2 2 1 1 2 2 ...
    %        1 2 1 2 ......
    
    combi=0*knots;
    
    % the easiest way to recover combi from knots is to proceed 1 row at a
    % time, i.e. one dimension at a time, and mark with a different label
    % (1,2,...P) all the equal points. We need of course as many labels 
    % as the number of different points in each dir!
    
    for dim=1:dim_tot
        
        % this is how many points per direction
        K=length(knots_per_dim{dim});
    
        % we start from a row of zeroes and we place 1....K in the right
        % poisitions by summations (each element of the row will be written
        % only once!)
        for k=1:K
            combi(dim,:) = combi(dim,:) + k*( knots(dim,:)==knots_per_dim{dim}(k) );
        end
        
    end
    
    
    % Now we can do the dot-multiplications among columns, the
    % multiplication by nodal values and the final sum! We proceed one
    % knot at a time
    
    for kk=1:S(i).size
        
        % dot-multiply all the lagrangian functions according to the
        % combi represented by the current knot
        
        f_loc=mono_lagr_eval{1}(:,combi(1,kk));
        for dim=2:dim_tot
            f_loc=f_loc.*mono_lagr_eval{dim}(:,combi(dim,kk));
        end

        % recover F, the corresponding value for the interpolating function in function_on_grid, with the global counter
        position = Sr.n(global_knot_counter);
        F_value = function_on_grid(position,:);
        
        % add the contribution of this knot to the sparse interpolation.
        % Its contribution is a matrix, since I am evaluating it on a bunch
        % of non_grid_points (one per row), and my function to be evaluated is vector-valued
        % which gives a bunch of columns. 
        
        f_values = f_values + coeff*f_loc*F_value;
        
        % in other words: f_values is the value of each output component on
        % each non_grid_point (hence matrix, as tall as non_grid_pts).
        %
        % we build it ``lagrange way'' by summing 
        %
        % (evaluations of each output component function in each grid point) x (lagrange polynomials evaluated in non_grid_points)
        % 
        % the lagrange polynomials are the same for every output component!
        % 
        % -) lagrange polynomials eval are in f_loc, it's a tall vector, as tall as non_grid_points. 
        % -) output components eval are in F_value, it's a fat vector, as fat as the number of outputs (it's a row of function_on_grid) 
        % 
        % so for each output component we want to compute
        %
        % f_loc * F_value(i) 
        %
        % is a tall vector, and we need to create the matrix [f_loc*F_value(1) f_loc*F_value(2) ....] i.e. vertical cat. 
        % we to this by [f_loc] * [F_value] (column vect by row vect multiplication)
        % 
        
        % update global counter
        global_knot_counter=global_knot_counter+1;
        
    end
    
    
end % of for loop on tensor grid





% % -------------------------------------------------
% % old code
% 
% % for each knot in knots I have to build the multidimensional lagrangian function and evaluate it
% 
% for idx_current_knot=1:S(i).size
%     
%     % recover F, the corresponding value for the interpolating function in function_on_grid, with the global counter
%     position = Sr.n(global_knot_counter);
%     F_value = function_on_grid(position);
%     
%     % this is the current knot where the lagrange function is centered
%     current_knot=knots(:,idx_current_knot);
%     
%     % compute the contribute of the current knot to the value of f in non_grid_points
%     f_values = f_values + coeff*F_value*lagr_eval_multidim(current_knot,knots_per_dim,non_grid_points);
%     
%     % update global counter
%     global_knot_counter=global_knot_counter+1;
% end
% 
% % old code
% % -------------------------------------------------

