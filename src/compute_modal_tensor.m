function U = compute_modal_tensor(S,S_values,domain,interval_map,flag)

% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,interval_map,'legendre')
% 
% or 
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,interval_map,'hermite')
%
% given a grid S in domain and the values on it, re-express it 
% as a modal expansion. 'flag' is optional, default is 'legendre' 
%
% output:
%
% U is a struct with fields U.size (the number of multi_indices), U.multi_indices, U.coeffs
%
%
% input:
%
% for LEGENDRE POLYNOMIALS
% 
% domain is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] 
% s.t. the interpolant is defined on (a1,b1) x (a2,b2) x ...  
%
% interval_map is a function that takes knots in S (the original grid, defined
% on a reference interval) to the knots in Sr (defined on domain). interval_map has to be defined like this:
%
%   interval_map =@(T) ....
%
% where input T is a matrix with one point per column (like the matrix [S.knots]) and the
% output is a matrix with the same size, and one point per column, shifted to domain
%
% for HERMITE POLYNOMIALS
%
% domain is a 2XN matrix = [mu1, mu2, mu3, ...; sigma1, sigma2, sigma3,...]
% the first variable has normal distribution with mean mu1 and std sigma1
% and so on...
% interval_map is a function that takes knots in S (the original grid,
% related to standard variables) to the knots in Sr (with mean and std in domain). interval_map has to be defined like this:
%
%   interval_map =@(T) ....
%
% where input T is a matrix with one point per column (like the matrix [S.knots]) and the
% output is a matrix with the same size, and one point per column, shifted
% to non standard distributions


if nargin==4
    flag='legendre';
end



% I will need the knots in each dimension separately. 
% As the number of knots is different in each direction, I use a cell array

nb_dim=size(S.knots,1);
knots_per_dim=cell(1,nb_dim);

for dim=1:nb_dim
    knots_per_dim{dim}=unique(S.knots(dim,:));
end

% The modal expansion in i-th direction uses up to degree k, with k as follows

degrees=zeros(1,nb_dim);
for dim=1:nb_dim
    degrees(dim) = length(knots_per_dim{dim})-1;
end

% the modal polynomials to be used are s.t. the corresponding multi-indices have
% all components less than or equal to the maximum degree

I = multiidx_box_set(degrees,0);

% return multiindex_set as 
U.multi_indices=I;

nb_multiindices=size(I,1);
U.size=nb_multiindices;

% safety check. I have solve a system, so I need the vandermonde matrix to be squared!
rows = S.size; % one equation per point
cols = nb_multiindices; % one unknown per polynomial

if rows~=cols, error('vandermonde matrix will not be square!'), end


% now create the vandermonde matrix with evaluation of each multi-index in every point
% of the grid

V = zeros(rows,cols);

for c=1:cols
    k = I(c,:);
    %vc = lege_eval_multidim(interval_map(S.knots),k,domain(1,:),domain(2,:));
    switch flag
        case 'legendre'
            vc = lege_eval_multidim(interval_map(S.knots),k,domain(1,:),domain(2,:));
        case 'hermite'
            vc = herm_eval_multidim(interval_map(S.knots),k,domain(1,:),domain(2,:));
        otherwise
            error('unknown family of polynomials')
    end    
    V(:,c)=vc';
end

% now solve the system

U.modal_coeffs = V \ S_values;



