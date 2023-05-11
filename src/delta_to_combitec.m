function I = delta_to_combitec(ii)


% DELTA_TO_COMBITEC computes the combination technique contributions of a multi-index
%
% I = delta_to_combitec(ii)
%
% takes as input a multi-index (row-vector) ii that describes a delta operator in the sparse grids construction
% and returns the set I of the corresponding quadrature/interpolation operators (lexicosorted). It does *not*
% return instead the coefficient of the index in the combination technique, whose value is trivial to compute, see below.
%
% In formulas,
%
% Delta^ii = \prod_{n=1}^N ( U_n^{ii_n} - U_n^{ii_n-1})
% 
% where U_n^{i_n} is the quadrature/interpolation operator along direction n at level i_n.
%
% The Delta^ii can be expressed as linear combination of tensorized U_n operators:
%
% Delta^ii = \sum_{kk \in K} c_kk \prod_n U_n^{kk_n}
%
% where c_kk = (-1)^sum(ii - kk);
%
% and I = sortrows(K),  while c_kk is not returned
% 
% For instance
%
% Delta^[2 3] = (U_1^2 -U_1^1) * (U_2^3 -U_2^2) = U_1^2 * U_2^3 - U_1^2 * U_2^2 - U_1^1 * U_2^3 + U_1^1 * U_2^2 
%
% and delta_to_combitec([2 3])
% 
% ans =
% 
%      1     2
%      1     3
%      2     2
%      2     3
% 
% Of course whenever ii_n = 1 for some n,  we should not take the difference ( U_n^{ii_n} - U_n^{ii_n-1}); delta_to_combitec does this too:
% 
% delta_to_combitec([2 2 3])
% 
% ans =
% 
%      1     1     2
%      1     1     3
%      1     2     2
%      1     2     3
%      2     1     2
%      2     1     3
%      2     2     2
%      2     2     3
%
% but delta_to_combitec([2 1 3])
% 
% ans =
% 
%      1     1     2
%      1     1     3
%      2     1     2
%      2     1     3

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


N = length(ii);

% I need to skip the components of ii equal to 1
N_eff = sum(ii>1); % other implementation have similar runtimes, like length(find(ii>1)) or sum(ii~=1) 

% preallocate I
I = ones(2^N_eff,N);

k_eff = 1;

% we work column-by-columns, i.e. loop over N
for k = 1:N
    
    % if the k-th component is not 1, replace the column of I, otherwise just leave the ones. The columns is obtained by repeating T times
    % the basic vector, which contains C copies of ii(k) and C copies of ii(k)-1.
    if ii(k) ~= 1
        
        times = 2^(k_eff - 1);
        copies  = 2^(N_eff - k_eff); % 2 x copies x times = #rows => 2 x copies x 2^(k_eff - 1)  = 2^N_eff => copies = 2^N_eff  / 2^k_eff
        
        row = 1;
        
        for it = 1 : times
            
            I(row : row + copies - 1,                                   k) = ii(k)-1; %<--- in this way, it's already lexicosrted. 
            I(                       row + copies : row + 2*copies - 1, k) = ii(k);
            row = row + 2*copies;
        end
        
        k_eff = k_eff + 1;
    end
    
end