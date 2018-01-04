function S = as_reduced(pts_list,wgs_list)

% S = AS_REDUCED(pts_list) returns a reduced sparse grids structure whose knots field is set equal to
% pts_list, while the other fields are []. Thus, isreduced(S)==true and
% S can then be used as input to e.g. PLOT_GRID or EVALUATE_ON_SPARSE_GRID. The knots in pts_list are supposed
% to be all different, i.e. the list is not checked for duplicates
%
% S = AS_REDUCED(pts_list,wgs_list) also set the weights field of S as wgs_list


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2017 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



S.knots = pts_list;
if nargin==2
    S.weights=wgs_list;
else
    S.weights=[];
end
S.n=[];
S.m=[];
