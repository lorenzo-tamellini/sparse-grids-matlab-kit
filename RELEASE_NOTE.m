%----------------------------------------------------------------------------------
% Sparse Grids Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------------------------------------

% RELEASE NOTES version 15.8
%
%
% -> 2015, Jun. 15  --> the new matlab toolbox syntax is now supported (see functions ACTIVATE_PARALLEL, CHECK_IF_PARALLEL_ON, CLOSE_PARALLEL)
%
%                   --> added MULTIIDX_GEN and MULTIIDX_BOX_SET to the tutorial
%
%
% -> 2015, Apr. 27  --> added function AS_REDUCED that converts a list of points / weights in a reduced grid
%                       structure
% 
%
% RELEASE NOTES version 14.12
%
%
% -> 2014, Dec. 2   --> CONVERT_TO_MODAL now handles the case of mixed polynomials in different directions. 
%
% -> 2014, Sep. 3   --> added support for Leja points (see KNOTS_LEJA, LEV2KNOTS_2STEP)
%
% -> 2014, Apr. 23  --> SMOLYAK_GRID, SMOLYAK_GRID_MULTIIDX_SET now return the multiidx associated to each tensor grid
% 
% -> 2014, Apr. 17  --> EVALUATE_ON_SPARSE_GRID now returns also the list of discarded points.
% 
%
%
% RELEASE NOTES version 14.4
%
%
% -> 2014, Apr. 09: --> LEV2KNOTS_NESTED has been renamed LEV2KNOTS_DOUBLING
%
%                   --> SMOLYAK_GRID_MULTIINDECES has been renamed SMOLYAK_GRID_MULTIIDX_SET. An error is thrown if the 
%                       multi-index set used is not sorted in lexicographic order
%
%                   --> SMOLYAK_GRID and SMOLYAK_GRID_MULTIIDX_SET have been rewritten and are now 2~3 times faster
%
%                   --> a new function, FAST_TD has been added to the package. It computes TD-like multiindex sets
%                       much faster than the standard MULTIIDX_SET
%
%                   --> it is now possible to recycle evaluations when using a sequence of sparse grids, see
%                       EVALUATE_ON_SPARSE_GRID, which supports Matlab parallel toolbox as well.
%                       The same features are also available for QUADRATURE_ON_SPARSE_GRID,  see help QUADRATURE_ON_SPARSE_GRID
%                       (inputs have been modified).
%
%                   --> it is now possible to mix random variables along different directions of a sparse grid,
%                       e.g. create a sparse grid with clenshaw-curtis points in one direction and gauss-hermite
%                       points in another, see help SMOLYAK_GRID, SMOLYAK_GRID_MULTIINDICES
%
%                   --> text ouput on screen is controlled by the boolean global variable MATLAB_SPARSE_KIT_VERBOSE
%                       (default value = 1 , i.e. output on screen is allowed)
%
%
%                   --> INTERPOLATE_ON_SPARSE_GRID, QUADRATURE_ON_SPARSE_GRID, CONVERT_TO_MODAL, COMPUTE_MODAL_TENSOR do not accept
%                       any longer MAP and WEIGHTS_COEFF input arguments, that should rather be passed as input to
%                       SMOLYAK_GRID, SMOLYAK_GRID_MULTIINDICES. 
%
%                   --> Input argument FLAG in CONVERT_TO_MODAL is now mandatory
%
%                   --> some of the inputs of INTERPOLATE_ON_SPARSE_GRID have been trasposed for consistency with the other
%                       functions. In particular, the size of FUNCTION_ON_GRID to be V x number_of_points_in_the_sparse_grid
%                       and consistently the size of the output matrix F_VALUES has been modified to V x (number_of_non_grid_points)
%                       with V such that f:R^N -> R^V. Also, the size of non_grid_points needs to be N x number_of_queried_evaluations 
%                       (that is, following the same convention as points in sparse grids)
%
%                   --> A warning is thrown if reduce_sparse_grid detects a tol which is inappropriately large 
%
%                   --> improved graphics handling for plot_grid
%
%                   --> help for many functions has been rewritten
%
%
% -> 2013, Oct. 24: INTERPOLATE_ON_SPARSE_GRID  interpolates vector-valued functions. 
%
% -> 2013, Jul. 03: added functions for CHEBYSHEV polynomials
%
% -> 2013, Apr. 18: fixed bug and backward compatibility for CHECK_SET_ADMISSIBILITY (becomes a standalone function on 2013, May 5)
% 
% -> 2012, Oct. 30: faster version of INTERPOLATE_ON_SPARSE_GRID
%
%
%
%
%
%
% sparse grid matlab kit august 2012 vs sparse grid toolkit november 2010
% 
% -----------------------------------------------
% 
% -> improved and expanded sparse grid tutorial
%
% -> new folder structure
%     /main includes the implementation of the functionalities of the code
%     /src includes inner code needed by functions in main 
%     /tools includes functions to generate knots, polynomials etc
%     /docs-examples contains tutorials and examples
% 
% -> the function check_set_admissibility has been added. It checks whether an index set satisfies the
%   admissibility condition
%
% -> /tools/polynomials now includes code to generate Chebyshev, Legendre, Lagrange and Hermite polynomials. 
%
% -> The functions to generate Lagrange polynomials are now called lagr_eval, lagr_eval_multidim instead
%    of monodim_lagrange_function, multidim_lagrange_function, 
%    to have the names of all the polynomial functions with the same structure
%     
% -> define_function_for_rule handles differently the case of isotropic rules, i.e. from
%       define_function_for_rule('TD',[1 1 1 1]) to define_function_for_rule('TD',4)
% 
% -> tools/knots_functions now includes code for KPN knots (nested knots for gaussian measure)
% 
% -> the functions to generate knots and weights for uniform and gaussian measures have been modified, from knots_lagrange/hermite/CC_st
%     to knots_uniform/gaussian/CC. It is now possible to specify the support interval of the rand vars, so that there is no need
%     to rescale knots and weights after having generated the grid. See the help for more details
%         
% -> however, for more complex cases (e.g. each rand var lives on its own interval), one still needs to
%    rescale afterwards. To this end, the function  get_interval_map.m      
%    will help generating the function that rescales the knots and the quadrature weights
% 
% -> knots_trap and midpoint have been removed
% 
% -> convert_to_modal.m function allows to compute the gPCE of a function starting from its sparse grid approximation.
%     The alogrithm implemented recasts the sum of lagrange polynomials into a sum of orthornormal polynomials and has 
%     been found to be more efficient than computing the gPCE coefficients by sparse grid quadrature.
%     Works for Chebyshev, Legendre and Hermite gPCE 
%
% -> multiindex_box_set.m becomes multiidx_box_set.m and the output order has been exchanged. Now the full box set is the first
%    output, instead of the second
