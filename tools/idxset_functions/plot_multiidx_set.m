function plot_multiidx_set(I,varargin)

% PLOT_MULTIIDX_SET(I) plots the index set I, in the case N=2 and N=3. 
%
% PLOT_MULTIIDX_SET(I,'PlotSpec',Spec_value, ...) specifies the plotting options to be used.

%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------



N = size(I,2);

switch N
    
    case 2
        if isempty(varargin)
            plot(I(:,1),I(:,2),'ok','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','k')
        else
            plot(I(:,1),I(:,2),varargin{:})
        end
        grid on
        maxI = max(max(I));
        axis([0 maxI+1 0 maxI+1])
        axis square
    
    case 3
        if isempty(varargin)
            plot3(I(:,1),I(:,2),I(:,3),'ok','LineWidth',2,'MarkerSize',12,'MarkerFaceColor','k')
        else
            plot(I(:,1),I(:,2),I(:,3),varargin{:})
        end
        grid on
        maxI = max(max(I));
        axis([0 maxI+1 0 maxI+1 0 maxI+1])
        axis square
        
    otherwise
        error('SparseGKit:WrongInput','cannot plot multidx with more than 3 dimensions')
end
