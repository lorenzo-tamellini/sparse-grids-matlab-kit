function L = piecewise_lin_eval_fast(central_knot,all_knots,non_grid_points,ngp_length)

% 

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2024 L. Tamellini, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------



nb_knots = length(all_knots);

switch nb_knots
    case 1
        L=ones(ngp_length,1);
        return
    case 2
        %error('piecewise function undefined for two knots')
    otherwise
        
end

L=zeros(ngp_length,1);

% idea of the code: assuming that non_grid_points is sorted, we only need to go through it once,
% and we check in which part of the hat function we are


%                                            /\
%                                           /  \ 
%                                          /    \
%                                         /      \
%                                        /        \
% x----------------x------------x-------p----o-----x----------x  <------ all_knots
%     *     *              *        *     * *|  *         * *
%     |                                  (i) |
%   non_grid_points                      central_knot, hat centered here
%
% |---------------- left part ---------------|---------------- right part ---------|



% evaluate the left-part (zero and increasing half ie _____/) of the hat function, but not if the hat is centered on the left border
% (central_knot = 1, hat shaped like this: \_____) 
if central_knot>1

    % preliminary check: if last of non_grid_points is smaller than p in the picture above, then we just return zeros
    if non_grid_points(end) < all_knots(central_knot - 1) % (@)
        return
    end

    % find first point in the support of the current hat function. L is zero until here, so no need to modify the
    % init value for previous entries. Note that the result of the search is never empty due to the check that we have done at line (@) 
    i = find(non_grid_points>=all_knots(central_knot - 1),1); 
    
    % From here, evaluate the first part of the support, until we are in it. Use a while to this end, increase i
    % one by one until non_grid_points(i) <= all_knots(central_knot) 
    
    % the slope of the hat
    m = 1/(all_knots(central_knot) - all_knots(central_knot-1));
    
    % the non_grid_points in the current part of the support. Careful: if central_knots = right border (hat shaped like this: ______/  ),
    % i.e. the current interval is at the end of the domain
    % the last point of non_grid_points is still in this interval, so we need to enforce that we do not ask for
    % non_grid_points(end+1)
    while i <= ngp_length && non_grid_points(i) <= all_knots(central_knot) 
        L(i) = m * (non_grid_points(i) - all_knots(central_knot-1));
        i = i+1;
    end
else
   % if central_knot = 1, we init the counter i to 1
   i = 1; 
end

% from here, evaluate the second part of the support (\___________), until we are in it. The same logic above applies. We need to
% make sure that:
% a) the hat function is not centered at the right border (central_knot = last knots, hat shaped like this: ______/ )
% b) if the hat is centered at the second-to-last knot (hat shaped like this: ______/\ ), i.e the current interval
%    is at the end of the domain, the last point of non_grid_points is still in this interval, so we need to enforce that we do not ask for
%    non_grid_points(end+1)

if central_knot<nb_knots
    
    m = 1/(all_knots(central_knot+1) - all_knots(central_knot));
    while i <= ngp_length && non_grid_points(i) < all_knots(central_knot+1)
        L(i) = m * ( - non_grid_points(i) + all_knots(central_knot+1));
        i = i +1;
    end
    
end

end