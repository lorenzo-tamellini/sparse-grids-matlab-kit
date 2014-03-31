function [A,i_ord]=mysortrows(A,Tol,i,i_ord,n)

% Similar to Matlab builtin function sortrows. Given a matrix A
% of real numbers, sorts the rows in lexicographic order; entries 
% that differ less than Tol are treated as equal (default Tol is 1e-14).
%
% usage: 
%   [B,i]=mysortrows(A,Tol)
%   input
%     A: input matrix
%     Tol: (optional) tolerance used to identify coincident entries 
%          (default 1e-14)
%   output
%     B: sorted matrix
%     i: index vector such that A(i,:)=B
%
% mysortrows has a recursive implementation.  
% recursive call: [A,i_ord]=mysortrows(A,Tol,i,i_ord,n)
% sorts the submatrix A(i,n:end); modifies the matrix itself and 
% stores the rows permutation in the global vector i_ord.



%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2014 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if nargin==1
    Tol=1e-14;
    n=1;
    i=1:size(A,1);
    i_ord=i;
elseif nargin==2
    n=1;
    i=1:size(A,1);
    i_ord=i;
end

%disp(['mysortrows level ',num2str(n)])

[tmp,ii]=sort(A(i,n));
A(i,:)=A(i(ii),:);
i_ord(i)=i_ord(i(ii));

if n<size(A,2)
    j=[1;diff(A(i,n))>Tol;1];
    jm=find(diff(j)==-1);
    jp=find(diff(j)==1);
    for k=1:length(jm)
        v=i(jm(k):jp(k))';
        [A,i_ord]=mysortrows(A,Tol,v,i_ord,n+1);
    end
end