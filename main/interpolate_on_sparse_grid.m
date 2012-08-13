function f_values = interpolate_on_sparse_grid(S,interval_map,Sr,function_on_grid,non_grid_points)

% value = interpolate_on_sparse_grid(S,interval_map,Sr,function_on_grid,non_grid_points)
%
% interpolates a function f defined on the sparse grid S. 
%
%   -> S is the original sparse grid
%   -> interval_map is a @-function that maps the interval containing knots of S to Sr:  Sr.knots=interval_map([S.knots]).
%           If no map is needed, set interval_map=[];
%   -> Sr is the reduced version of S
%   -> function_on_grid is a vector containing the evaluation of the function on the points of the (reduced) grid
%   -> non_grid_points is the set of points (not belonging to the sparse grid) where I want to evaluate the function.
%           non_grid_points is a matrix, each row is a different point


f_values=zeros( size(non_grid_points,1) , 1);

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
    
    % for each knot in knots I have to build the multidimensional lagrangian function and evaluate it
    
    for idx_current_knot=1:S(i).size
        
        % recover F, the corresponding value for the interpolating function in function_on_grid, with the global counter 
        position = Sr.n(global_knot_counter);
        F_value = function_on_grid(position);

        % this is the current knot where the lagrange function is centered
        current_knot=knots(:,idx_current_knot);
        
        % compute the contribute of the current knot to the value of f in non_grid_points 
        f_values = f_values + coeff*F_value*lagr_eval_multidim(current_knot,knots_per_dim,non_grid_points);
        
        % update global counter
        global_knot_counter=global_knot_counter+1;
    end
    
end % of for loop on tensor grid

end
