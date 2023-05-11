function h= plot3_sparse_grid(S,dims,varargin)

% h = plot3_sparse_grid(S,dims,varargin)
%
% PLOT3_SPARSE_GRID(S) plots S, which is a sparse grid in 3D. S can be either reduced or not. S can also be a tensor grid
%       (use is_tensor_grid(S) to verify, type HELP IS_TENSOR_GRID for details)
%
% PLOT3_SPARSE_GRID(S,[d1 d2 d3]) plots the components d1, d2, d3 of the points in S if S is more than 3D
%
% PLOT3_SPARSE_GRID(S,[d1,d2,d3],varargin) defines the style of the plot, e.g.
%
% plot3_sparse_grid(S,[],'color','g','marker','o')
%
% h=PLOT3_SPARSE_GRID(S)
%
% returns the handle to the plot


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



if ~exist('dims','var') || isempty(dims) 
    dims=[1 2 3];
end

x=[S.knots];

if nargin==1 || nargin==2
    % use a default plot style
    h=plot3(x(dims(1),:),x(dims(2),:),x(dims(3),:),'ok','MarkerFaceColor','k');
else
    % use style provided
    h=plot3(x(dims(1),:),x(dims(2),:),x(dims(3),:),varargin{:},'LineStyle','none');
end
grid on


