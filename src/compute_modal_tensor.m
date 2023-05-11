function U = compute_modal_tensor(S,S_values,domain,flags)

% COMPUTE_MODAL_TENSOR given a tensor grid and the values on it, re-express the interpolant 
% as a modal expansion.
% 
% U=COMPUTE_MODAL_TENSOR(S,S_VALUES,DOMAIN,'legendre') considers the tensor grid S on the
%       hyper-rectangle DOMAIN with associated point evaluations S_VALUES and converts the
%       resulting lagrangian multivariate interpolant to a sum of Legendre polynomials.
%       DOMAIN is a 2xN matrix = [a1, a2, a3, ...; b1, b2, b3, ...] 
%       defining the hyper-rectangluar domain of the sparse grid: (a1,b1) x (a2,b2) x ...  
%       U is a struct with fields U.size (the number of Legendre polynomials needed), 
%       U.multi_indices (one multi-index per Legendre polynomial), U.coeffs
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'chebyshev') works as the previous call, using 
%       Chebyshev polynomials
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'hermite') works as the previous call, using 
%       Hermite polynomials. Here DOMAIN is a 2XN matrix = [mu1, mu2, mu3, ...; sigma1, sigma2, sigma3,...]
%       such that the first variable has normal distribution with mean mu1 and std sigma1
%       and so on.
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'laguerre') works as the previous call, using 
%       Laguerre polynomials. Here DOMAIN is a 1XN matrix = [lambda1, lambda2, lambda3, ...]
%       such that the first variable has exponential distribution with parameter lambda1 and so on.
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'generalized laguerre') works as the previous call, using 
%       generalized Laguerre polynomials. Here DOMAIN is a 2XN matrix = [alpha1, alpha2, alpha3, ...; beta1, beta2, beta3, ...]
%       such that the first variable has Gamma distribution with parameters alpha1 and beta1 and so on.
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,'jacobi') works as the previous call, using 
%       Jacobi polynomials. Here DOMAIN is a 4XN matrix = [alpha1, alpha2, alpha3, ...; beta1, beta2, beta3, ...; a1, a2, a3, ...; b1, b2, b3, ...]
%       such that the first variable has Beta distribution with parameters alpha1 and beta1 on the interval [a1,b1] and so on.
%
% U=COMPUTE_MODAL_TENSOR(S,S_values,domain,{<family1>,<family2>,<family3>,...}) works as the previous call, using 
%       polynomials of type <family-n> in direction n. Here DOMAIN is a cell array of lenght N 
%       where each cell contains the matrix giving the parameters of the n-th family of polynomials. 
%       For example:
%
%       COMPUTE_MODAL_TENSOR(S,S_VALUES,DOMAIN,{'legendre','hermite','laguerre','jacobi','legendre'})
%       with 
%       DOMAIN = {[a1;b1], [mu1;sigma1], lambda, [alpha2;beta2;a2;b2], [a3;b3]}

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if any(~ismember(flags,{'legendre','chebyshev','hermite','laguerre','generalized laguerre','jacobi'}))
    error('SparseGKit:WrongInput',['Input argument FLAG unrecognized. '...
        'Please note that COMPUTE_MODAL_TENSOR does not accept INTERVAL_MAP '...
        'input argument any longer. '...
        'Type help convert_to_modal for help. '...
        'This error message will not be shown in future releases of SPARSE-GRID-MATLAB-KIT']);
end

if iscell(flags) && ~iscell(domain) 
    errmsg = ['Input argument DOMAIN must be a cell array. ' ... 
            'Please note that CONVERT_TO_MODAL has been changed after release 18.10. ' ...
            'The domain for the case of polynomials of "mixed" type is now a cell array, each cell containing the domain for the corresponding polynomial. ' ...
            'Type help convert_to_modal for help. '...
            'This message will disappear from future relesases of SPARSE-GRID-MATLAB-KIT.'];
    error('SparseGKit:WrongInput',strcat(errmsg));
end

% I will need the knots in each dimension separately. 
% As the number of knots is different in each direction, I use a cell array

nb_dim=size(S.knots,1);

% knots_per_dim=cell(1,nb_dim);
% for dim=1:nb_dim
%     knots_per_dim{dim}=unique(S.knots(dim,:));
% end

knots_per_dim=S.knots_per_dim;


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

if rows~=cols, error('SparseGKit:FailedSanityChk','vandermonde matrix will not be square!'), end


% now create the vandermonde matrix with evaluation of each multi-index in every point
% of the grid

V = zeros(rows,cols);

if length(flags)==1 || ischar(flags) % the second condition for when the function is called on one sigle family of polynomials
    for c=1:cols
        k = I(c,:);
        %vc = lege_eval_multidim(interval_map(S.knots),k,domain(1,:),domain(2,:));
        switch flags
            case 'legendre'
                vc = lege_eval_multidim(S.knots,k,domain(1,:),domain(2,:));
            case 'hermite'
                vc = herm_eval_multidim(S.knots,k,domain(1,:),domain(2,:));
            case 'chebyshev'
                vc = cheb_eval_multidim(S.knots,k,domain(1,:),domain(2,:));
            case 'laguerre'
                vc = lagu_eval_multidim(S.knots,k,domain);
            case 'generalized laguerre'
                vc = generalized_lagu_eval_multidim(S.knots,k,domain(1,:),domain(2,:));
            case 'jacobi'
                vc = jacobi_prob_eval_multidim(S.knots,k,domain(1,:),domain(2,:),domain(3,:),domain(4,:));
            otherwise
                error('SparseGKit:WrongInput','unknown family of polynomials')
        end
        V(:,c)=vc';
    end
else
    for c=1:cols
        k = I(c,:);
        vc=ones(1,size(S.knots,2));
        for n=1:nb_dim
            switch flags{n}
                case 'legendre'
                    % vc = vc.*lege_eval(S.knots(n,:),k(n),domain(1,n),domain(2,n));
                    vc = vc.*lege_eval(S.knots(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'hermite'
                    % vc = vc.*herm_eval(S.knots(n,:),k(n),domain(1,n),domain(2,n));
                    vc = vc.*herm_eval(S.knots(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'chebyshev'
                    % vc = vc.*cheb_eval(S.knots(n,:),k(n),domain(1,n),domain(2,n));
                    vc = vc.*cheb_eval(S.knots(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'laguerre'
                    vc = vc.*lagu_eval_multidim(S.knots(n,:),k(n),domain{n});
                case 'generalized laguerre'
                    vc = vc.*generalized_lagu_eval_multidim(S.knots(n,:),k(n),domain{n}(1),domain{n}(2));
                case 'jacobi'
                    vc = jacobi_prob_eval_multidim(S.knots(n,:),k(n),domain{n}(1),domain{n}(2),domain{n}(3),domain{n}(4));
                otherwise
                    error('SparseGKit:WrongInput','unknown family of polynomials')
            end
        end
        V(:,c)=vc';
        
    end
end

% now solve the system

U.modal_coeffs = V \ S_values;
