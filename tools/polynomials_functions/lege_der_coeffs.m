function [dc,dercoeffs] = lege_der_coeffs(degmax,a,b)


% LEGE_DER_COEFFS returns the coefficients of the derivatives of orhtonormal legendre polynomials, 
% ordered according matlab convention for polynomials (i.e. coefficient of degree 0 at the rightmost place)
%
% DC = LEGE_DER_COEFFS(DEG,A,B) returns the coefficients of the derivative of the legendre polynomial 
%       of degree DEG orthonormal with respect to the weight 1/(b-a)
%
% [DC, DERCOEFFS] = LEGE_DER_COEFFS(DEG,A,B) returns the coefficients of the derivative of **all** the polynomials
%       with degree <= DEG, stored by row in a matrix POLCOEFFS.
%       The i-th row of POLCOEFFS contains the coefficients of the (i-1)-degree polynomial,
%       i.e. the constant polynomial in the first row, linear polynomial in the second row etc


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2017 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if degmax==0
    dc=[];
    dercoeffs=[];
    return
end

% compute coefficients of polynomials
[~,polcoeffs] = lege_coeffs(degmax,a,b);

%initializations
nbpol=size(polcoeffs,1);
degmaxcheck=size(polcoeffs,2)-1;

if degmaxcheck~=degmax
    error('degmaxcheck~=degmax')
end

dercoeffs=zeros(nbpol,degmax+1);

% for each polynomial, derive using built-in matlab functions and store it
% in dercoeffs
for p=1:nbpol
    dx = polyder(polcoeffs(p,:));
    dx_size = length(dx);
    dercoeffs(p,end-dx_size+1:end) = dx;
end

dc = dercoeffs(end,:);