function [lc,polcoeffs] = lege_coeffs(degmax,a,b)

% LEGE_COEFFS returns the coefficients of orhtonormal legendre polynomials, 
% ordered according matlab convention for polynomials (i.e. coefficient of degree 0 at the rightmost place)
% 
% LC = LEGE_COEFFS(DEG,A,B) returns the coefficients of the univariate polynomial of degree DEG
%       orthonormal with respect to the weight 1/(B-A), contained in the vector LC
% 
% [LC, POLCOEFFS] = LEGE_COEFFS(DEG,A,B) returns the coefficients of **all** the polynomials
%       with degree <= DEG, stored by row in a matrix POLCOEFFS.
%       The i-th row of POLCOEFFS contains the coefficients of the (i-1)-degree polynomial,
%       i.e. the constant polynomial in the first row, linear polynomial in the second row etc

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2017 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



% we begin with a temp matrix where the i-th row has the coefficients of the (i-1)-th
% polynomial, ordered in the opposite order than Matlab, i.e. 
% from left to right FROM MIN TO MAX (easier to code, flip at the end). 
tpolcoeffs=zeros(degmax+1,degmax+1);

% init with L0=1, L1=x
tpolcoeffs(1,1)=1;
tpolcoeffs(2,2)=1; 

% recursive step for standard Legendre polynomials in -1,1:
%
%  L_{n} = (2*n - 1) / n * x .* L_{n_1} - (n-1) / n * L_{n-2};
%
% note that here subindex n denotes the degree of the polynomial, n
% starting from 0
%  
%
% for loop 
for kap=2:degmax % kap is the polynomial degree, same as n above
    
    %pos is the shift due to the numbering of rows of matlab
    pos=kap+1; 
    
    % the first term of the recursion multiplies the (n-1)-th polynomial by
    % x, so all of its coefficients now refers to a monomial with 1 degree
    % more. In practice, we need to shift their position 1 place rightward
    % on the row. 
    tpolcoeffs(pos,:)= [ 0  (2*kap - 1)/kap*tpolcoeffs(pos-1,1:pos-1)  zeros(1,degmax-kap) ]; 
    
    % the second term involves the n-2-th polynomial
    tpolcoeffs(pos,:) = tpolcoeffs(pos,:) -  (kap - 1)/kap*tpolcoeffs(pos-2,:);
end

% now orthonormalize, by multiplying P_n * sqrt(2*n+1)
for kap=0:degmax 
    pos=kap+1; 
    tpolcoeffs(pos,:) = tpolcoeffs(pos,:) * sqrt(2*kap+1); 
end


% now flip to get back to matlab notation
tpolcoeffs=fliplr(tpolcoeffs);


% take into account different intervals. 
% we have to accomodate each polynomial with the change of variable, one
% monomial at a time. In particular, for each monomial c*x^m in each
% polynomial, we have to generate the polynomial corresponding to
% c* ( 2/(b-a) x + (-b-a)/(b-a) )^m, and add this contribution to the final 
% value of the polynomial. We can precompute these power polynomials
% once for all before the computation starts. Each has a different length
% so we store them in a cell-array


c1=2/(b-a);
c2=-(b+a)/(b-a);

pws = cell(1,degmax);
pw=[c1 c2];

pws{1}=pw;
for dd=2:degmax %starting from square power
    pw=conv(pw,[c1,c2]); % this computes the product of the polynomial in conv with the polynomial represented by [c1 c2]
    pws{dd}=pw;
end



% now apply the transformation. the 0-th polynomial (constant) does not change 
polcoeffs=zeros(size(tpolcoeffs));
polcoeffs(1,end)=tpolcoeffs(1,end);


% scroll the rest of the polynomials
for kap=1:degmax 

    pos=kap+1;

    % the constant term of the kap-th polynomial does not change anyway
    polcoeffs(pos,end)=tpolcoeffs(pos,end);


    % now consider the monomials one at a time and apply the change of
    % variables
    
    for mon=1:kap % the degree of the kap-th polynomial are kap at most
    
        % now we are working in the Matlab reference system, so that the
        % monomials are ordered from right to left
        monpos = degmax + 1 - mon ;

        % recover the change of variable contribution, multiply by the
        % initial coefficient and add zeros in front of pw if needed, 
        % so that we can then sum it to polcoeff_fin    
        pw=pws{mon};
        pw_len = length(pw);
        pw = tpolcoeffs(pos,monpos) * [zeros(1,degmax+1-pw_len) pw];
        
        polcoeffs(pos,:)= polcoeffs(pos,:) + pw;

    end

end


% assign output
if degmax>0
    lc=polcoeffs(end,:);
else
    lc=polcoeffs(1,:); % this is because we have 2 polynomials in any case, given that we have 2 polynomials in the recursion
end