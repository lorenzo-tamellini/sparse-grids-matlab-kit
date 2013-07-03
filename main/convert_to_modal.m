function [modal_coeffs,K] = convert_to_modal(S,Sr,nodal_values,domain,interval_map,flag)

% [modal_coeffs,K] = CONVERT_TO_MODAL(S,Sr,nodal_values,domain,interval_map)
%
% recast a sparse grid interpolant as a sum of orthogonal polynomials
% i.e. computes the spectral expansion of the interpolant.
%
% the multi-indices of the polynomials appearing in the conversion are
% stored in matrix K, one row per multi-index
%
% flag='legendre','hermite','chebyshev' is optional. Default is legendre
%
% -> S is the original grid
% -> Sr is the reduced grid
% -> nodal_values are the values of the function on the reduced sparse grid
%
% for LEGENDRE and CHEBYSHEV polynomials
%
% -> domain is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] 
%   s.t. the interpolant is defined on (a1,b1) x (a2,b2) x ...  
%
% -> interval_map is a function that takes knots in S (the original grid, defined
%   on a reference interval) to the knots in Sr (defined on domain). interval_map has to be defined like this:
%
%   interval_map =@(T) ....
%
%   where input T is a matrix with one point per column (like the matrix [S.knots]) and the
%   output is a matrix with the same size, and one point per column, shifted to domain. If no map is
%   needed, set interval_map=[];
%
%
% for HERMITE polynomials
%
% -> domain is a 2XN matrix = [mu1, mu2, mu3, ...; sigma1, sigma2, sigma3,...]
%   the first variable has normal distribution with mean mu1 and std sigma1
%   and so on...
%
% -> interval_map is a function that takes knots in S (the original grid,
%   related to standard variables) to the knots in Sr (with mean and std in domain). interval_map has to be defined like this:
%
%   interval_map =@(T) ....
%
%   where input T is a matrix with one point per column (like the matrix [S.knots]) and the
%   output is a matrix with the same size, and one point per column, shifted
%   to non standard distributions. If no map is needed, set interval_map=[];


% fix input 
if nargin==5
    flag='legendre';
end
if isempty(interval_map)
    interval_map = @(t) t;
end


% nodal values has to be column
[r,c]=size(nodal_values);
if r<c
    nodal_values=nodal_values';
end

% N is be the number of random variables
N=size(Sr.knots,1);


% one tensor grid at a time, compute the modal equivalent, then sum up
% how many grid to process?
nb_tensor_grids = size(S,2);


% each tensor grid will result in vector of coefficients, u, that
% are the coefficients of the modal (Legendre,Hermite) expansion of the 
% lagrangian interpolant. To combine these coefficients properly, I need 
% to store not only the coefficients u themselves, but also the associated
% multi-indices. I store this information into a structure array, U, 
% whose lenght is nb_tensor_grids, with fields U(i).modal_coeffs, U(i).multi_indices, U(i).size (i.e. 
% how many polynomials are associated to each grid)

% I will need this counter
global_values_counter=0; 

% I also do a trick to preallocate U
U(nb_tensor_grids).multi_indices=[];
U(nb_tensor_grids).size=[];
U(nb_tensor_grids).modal_coeffs=[];

for g = 1: nb_tensor_grids
    
    % some of the grids in S structure are empty, so I skip them
    if isempty(S(g).knots)
        continue
    end

    % recover the values of f on the current grid. 
    % To do this, I need to use Sr.n, that is the mapping from [S.knots] to Sr.knots
    % To locate the knots of the g-th grid in the mapping Sr.n, i use global_values_counter    
    grid_knots = global_values_counter + 1 : ( global_values_counter + S(g).size );
    
    % extract the values of the g-th grid
    S_values = nodal_values(Sr.n(grid_knots));

    % then update global_values_counter
    global_values_counter = global_values_counter + S(g).size ;
        
    % and compute modal
    U(g) = compute_modal_tensor(S(g),S_values,domain,interval_map,flag);
    
end

% now I need to put together all the U.multi_indices, summing the coefficients that
% share the same multi_index. Same procedure as reduce sparse grid

% First I build a list of all multi_indices and coefficients
% I preallocate them to speed up code

tot_midx=sum([U.size]);

% this is the container of multi_indices
All = zeros(tot_midx,N);

% this is the container of coefficients
all_coeffs = zeros(tot_midx,1);

% this is the index that scrolls them
l=0; 

% fill All and all_coeffs
for g = 1: nb_tensor_grids

    if isempty(S(g).knots)
        continue
    end

    All(l+1:l+U(g).size,:) = U(g).multi_indices;
    all_coeffs(l+1:l+U(g).size) = U(g).modal_coeffs*S(g).coeff;
    l=l+U(g).size;
end


% get rid of identical ones, by sorting and taking differences between cosecutive rows
% I need to store the sorter vector, to sum up all coefficients with
% the same multi_index

[All_sorted,sorter] = sortrows(All);

% taking differences between consecutive rows, I have a way
% to build the unique version of All. Remember diff has 1 row
% less, so I have to add a ficticious row at the end of All

dAll_sorted=diff( [All_sorted; ...
                   All_sorted(end,:)+1] );

% if a row of dAll_sorted is made of 0s then I have
% to discard the corresponding row of All

selector = sum(abs(dAll_sorted),2); % use abs(dAll_sorted), difference could be negative!
selector = selector>0;

% store uniqued version of All
K = All_sorted(selector,:);

% now I have to sum all the coefficients.

% sort coefficients in the same way as multi_indices
all_coeffs_sorted = all_coeffs(sorter);

% initialize modal_coeffs, that is the final output
nb_modal_coeffs = size(K,1);
modal_coeffs = zeros(nb_modal_coeffs,1);

% compute them one at a time

l=1; %scrolls all_coeffs_sorted
for k= 1 : nb_modal_coeffs
    ss=all_coeffs_sorted(l);
    while selector(l)~=1
        l=l+1;
        ss = ss+all_coeffs_sorted(l);
    end
    modal_coeffs(k) = ss; 
    l=l+1;
end