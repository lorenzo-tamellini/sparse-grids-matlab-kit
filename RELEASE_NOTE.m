%----------------------------------------------------------------------------------
% Sparse Grids Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------------------------------------
% 
%
%
% RELEASE NOTES version 23.05 (Robert)
%
%
% -> 2023, May 7   --> added testing unit in docs-examples/testing_unit
% 
% -> 2023, May 7   --> added function KNOTS_TRIANGULAR_LEJA i.e., weighted Leja for quadrature with respect to the
%                      triangular pdf
% 
% -> 2023, May 6   --> renamed SMOLYAK_GRID, SMOLYAK_GRID_MULTIIDX_SET, SMOLYAK_GRID_ADD_MULTIIDX, SMOLYAK_GRID_QUICK_PRESET to CREATE_SPARSE_GRID,
%                      CREATE_SPARSE_GRID_MULTIIDX_SET, CREATE_SPARSE_GRID_ADD_MULTIIDX, CREATE_SPARSE_GRID_QUICK_PRESET 
% 
%                  --> renamed ISSMOLYAK to IS_SPARSE_GRID to avoid using misleading function names and ISTENSOR to IS_TENSOR_GRID 
%                      for consistency
%
%
% RELEASE NOTES version 22.02 (California)
%
%
% -> 2022, Feb. 16 --> added function SMOLYAK_GRID_QUICK_PRESET, to generate in one command a simple sparse grid
%                      (CC knots in [-,1,1], lev2knots_doubling, idxset: @(ii) sum(ii-1) \leq w) and its reduced version
%
% -> 2022, Feb. 14 --> added PLOT3_SPARSE_GRID to plot 3-dimensional sparse grids (or 3-dimensional cuts of sparse grids with N>3)
%
% -> 2022, Feb. 11 --> added UNIVARIATE_INTERPOLANT to compute the one-dimensional Lagrangian interpolant of a scalar-valued function.
%                   
%                  --> added TENSOR_TO_SPARSE that converts a tensor grid structure into a sparse grid structure.
%                      This is useful when computing quadrature and interpolation on tensor grids using 
%                      the function evaluate/quadrature_on_sparse_grid, that contain a check whether the input is
%                      a sparse grid or not.
%                 
% -> 2022, Feb. 1   --> KNOTS_GK returns now GK points and weights for integration 
%                       w.r.t. general Gaussian densities with mean mi and standard deviation sigma.
%
%                   --> compatibility with GNU Octave have been partially tested. 
%
% -> 2021, Dec 1    --> added function KNOTS_BETA_LEJA, i.e. weighted Leja points (and their symmetric version) for quadrature with respect to the 
%                       Beta weight function with range [x_a,x_b], and parameters alpha,beta>-1. 
%                       Added example file where Beta Leja quadrature and interpolation convergence tests for different univariate 
%                       knots for Beta random variables are compared, TEST_CONVERGENCE_BETA_LEJA.m
%
%                   --> added function KNOTS_GAMMA_LEJA, i.e. weighted Leja points for quadrature with respect to the 
%                       Gamma weight function with parameter alpha>-1, beta=1. 
%                       Added example file where Gamma Leja quadrature and interpolation convergence tests for different univariate 
%                       knots for Gamma random variables are compared, TEST_CONVERGENCE_GAMMA_LEJA.m
%
%                   --> added function KNOTS_EXPONENTIAL_LEJA, i.e., weighted Leja for quadrature with respect to the
%                       exponential weight function. Added example file where exponential Leja are computed
%                       and then quadrature and interpolation convergence tests for different univariate 
%                       knots for exponential random variables are compared, TEST_COMPUTE_EXPONENTIAL_LEJA_AND_CONVERGENGE_TEST.m
%
%                   --> extend CONVERT_TO_MODAL and COMPUTE_MODAL_TENSOR to Laguerre, generalized Laguerre and Jacobi polynomials. 
%                       Note that the input DOMAIN for the case of polynomials of "mixed" type is now a cell array, 
%                       each cell containing the domain for the corresponding polynomial. 
%                       Corresponding tests have been added to TEST_CONVERT_TO MODAL.   
%
%                   --> added function LAGU_EVAL and LAGU_EVAL_MULTIDIM to generate one-dimensional and multi-dimensional Laguerre polynomials 
% 
%                   --> added function GENERALIZED_LAGU_EVAL and GENERALIZED_LAGU_EVAL_MULTIDIM to generate one-dimensional and multi-dimensional generalized Laguerre polynomials
%
%                   --> added function JACOBI_EVAL and JACOBI_EVAL_MULTIDIM to generate one-dimensional and multi-dimensional Jacobi polynomials
%                   
%                   --> added function KNOTS_BETA to generate Gauss-Jacobi knots and weights for integration
%                       w.r.t. Beta distributions.
%                   
%                   --> added function KNOTS_GAMMA to generate Gauss-generalized Laguerre knots and weights for integration
%                       w.r.t. Gamma distributions.
% 
%                   --> added function KNOTS_EXPONENTIAL to generate Gauss-Laguerre knots and weights for integration
%                       w.r.t. exponential distributions.
%
%                   --> KNOTS_GAUSSIAN_LEJA has been renamed to KNOTS_NORMAL_LEJA and now returns both standard and symmetric Leja points+weights for integration 
%                       w.r.t. general Gaussian densities with mean mi and standard deviation sigma. 
%                       The function requires now four inputs: n (the number of points required), mi, and sigma of
%                       the pdf, and the type of knots ('line'/'sym_line')
%
%                   --> KNOTS_GAUSSIAN has been renamed to KNOTS_NORMAL, for naming consistency
%
% -> 2020, Jul. 21  --> for robustness, SMOLYAK_GRID, SMOLYAK_GRID_MULTIIDX_SET no longer accepts MAP and COEFF_WEIGHTS as inputs. The knots and weights
%                       need to be already properly rescaled beforehand, by using the optional arguments to defind correctly the 1D families of points to be used, 
%                       and then invoking e.g.
%                       [S,C] = SMOLYAK_GRID(N,W,{@knots1, @knots2, ...},{@m1, @m2 ...}),
%
%                   --> removed QUADRATURE_ON_SPARSE_GRID_LEGACY and EVALUATE_ON_SPARSE_GRIDS_LEGACY
%
%                   --> improved the help function of PLOT_SPARSE_GRIDS_INTERPOLANT and COMPUTE_SOBOL_INDICES_FROM_SPARSE_GRID, DERIVE_SPARSE_GRID.
%
%                   --> added function PLOT_MULTIIDX_SET. Works for N=2 and N=3 only. For larger dimensions, the user needs to specify which dimensions should be plotted
%
%                   --> PLOT_SPARSE_GRID_INTERPOLANT does not generate a new figure, unless cuts are requested  (N>3)
%
%                   --> for convenience, QUADRATURE_ON_SPARSE_GRID now can also be called as 
%
%                       QUADRATURE_ON_SPARSE_GRID(F_VALS,SR)                   
%   
%                       where F_VALS is a row vector containing the values of the function whose quadrature we need to compute, that we might have already available. 
%                       It works also for multivariate functions (each component as row of a matrix)
%
%                   --> added ISEQUAL_TOL, ISLEXICO_TOL, FIND_LEXICOGRAPHIC_TOL, to deal with vectors with numerical noise
%
%                   --> added HESSIAN_SPARSE_GRID to compute the hessian of a function by taking finite differences of its sparse grids approximation
%
%                   --> ADAPT_SPARSE_GRID can now work with different families of knots in different directions. However, there are some restrictions, and
%                       the user is encouraged to read the help function and take a look at the new version of TUTORIAL_ADAPTIVE.M
%
%                   --> clarified that the option "buffer" for ADAPT_SPARSE_GRID can only be used in certain cases.
%                       The user is encouraged to read the help function and take a look at the new version of TUTORIAL_ADAPTIVE.M 
%                       The default value of "controls.var_buffer_size" has been changed to min(N_full,5) to N_full.
% 
%                   --> fixed a bug in ADAPT_SPARSE_GRID that would halt with error for certain values of controls.pts_tol   
%
%                   --> for consistency with the literature,  the knots KPN (Kronrod Patterson Normal, nested knots for gaussian rand var) have been renamed
%                       as Genz--Keister. In particular: 
%                       -- KNOTS_KPN is now KNOTS_GK; 
%                       -- LEV2KNOTS_KPN is now LEV2KNOTS_GK; 
%                       -- KPN_LEV_TABLE is now GK_LEV_TABLE 
%               
%                   --> renamed DETECT_UNSUFFICIENT_TOLERANCE to DETECT_INSUFFICIENT_TOLERANCE
%                       
% 
% -> 2020, Apr. 1   --> fixed a bug in SMOLYAK_GRID_MULTIIDX_SET, which would otherwise throw an error when called as SMOLYAK_GRID_MULTIIDX_SET(C,KNOTS,LEV2KNOTS,[])
%
%                   --> Added more NORMAL_LEJA (formerly GAUSSIAN_LEJA) points, we now have 150 instead of 50. Also, they are no more saved as a .mat file, 
%                       but pasted in ascii format in the knots_gaussian_leja, it's much faster than loading the matlab file each time
%
%
% -> 2020, Mar. 22  --> faster version of ADAPT_SPARSE_GRID, by using the new function SMOLYAK_GRID_ADD_MULTIIDX (see below) 
%
%                   --> added functions to test whether two tensor grids and two sparse grids are equal, ISEQUAL_SPARSE_GRIDS and ISEQUAL_TENSOR_GRIDS
%
%                   --> added functions to add a single multi-idx to a sparse grid, SMOLYAK_GRID_ADD_MULTIIDX and the ancillary function DELTA_TO_COMBITEC	
%
%                   --> added function to compute the combination technique coefficients from a multiidx set, COMBINATION TECHNIQUE	
%
%                   --> changed interface of function ISTENSOR. The output can now take 3 values: 1 if the input is a tensor grid, -1 if the input is
%                       a tensor grid stored in a sparse grid struct (i.e., a standard tensor with the additional fields idx and coeff) and 0 otherwise	
%
%
% -> 2019, Feb. 23  --> added midpoint and trapezoidal univariate quadrature/intepolation rules, and LEV2KNOTS_TRIPLING 
%                       to obtain nested sequences of midpoint knots
%
%
% -> 2019, Feb. 7   --> fixed a low-level bug in evaluate on sparse grid that would make the function stop with an error (very rare)
%
%
% -> 2019, Feb. 7   --> added error messages for functions that have been renamed
%
%
% -> 2019, Feb. 7   --> the output of ADAPT_SPARSE_GRID now contains a field "nested" which is set to TRUE if nested points were used
%
%
%
% RELEASE NOTES version 18.10 (Esperanza)
%
%
% -> 2018, Oct. 10  --> added DERIVE_SPARSE_GRID to compute gradients of a sparse grid interpolant
%
%                   --> added COMPUTE_SOBOL_INDICES_FROM_SPARSE_GRID to compute Sobol indices of a function by its sparse grid interpolant
%
%                   --> added KNOTS_GAUSSIAN_LEJA, i.e., weighted Leja for quadrature with respect to the
%                       gaussian weight function. Added example file where Gaussian Leja are computed
%                       and then quadrature and interpolation convergence tests for different univariate 
%                       knots for Gaussian random variables are compared, TEST_COMPUTE_GAUSSIAN_LEJA_AND_CONVERGENGE_TEST.m
%
%                   --> added EXPORT_SPARSE_GRID_TO_FILE to save a reduced grid on file (knots and weights)
%
%                   --> renamed PLOT_GRID to PLOT_SPARSE_GRID for naming consistency
%
%                   --> added PLOT_SPARSE_GRID_INTERPOLANT to plot the sparse grid approximation of a function.
%                       (different plots for the cases N=2, N=3,  N>3)
% 
%                   --> added summary at the beginning of tutorial
%
%                   --> added slides about the software as additional documentation
%
%                   --> added field "size" to reduced sparse grid, that contains the number of different points in the sparse grid,
%                       i.e.,  Sr.size == length(Sr.weights) == size(Sr.knots,2)
%                       
%                   --> renamed AS_REDUCED -> ASREDUCED, for consistency with other names
%
%                   --> fixed bug in COMPARE_SPARSE_GRIDS with use of nargin
%
%                   --> fixed bug in PLOT_SPARSE_GRID for nargin == 2
%
% -> 2018, Jun. 28  --> added example of convergence study for interpolation error with adaptive sparse grids,
%                       see tutorial_adaptive.m
%
%
%
% RELEASE NOTES version 17.5 (Trent)
%
%
% -> 2017, May  17  --> added ADAPT_SPARSE_GRIDS and TUTORIAL_ADAPTIVE and PLOT_IDX_STATUS to show how adapt
%                       sparse explores 2-dimensional spaces
%
% -> 2017, Feb. 14  --> CONVERT_TO_MODAL now can operate on vector-valued functions. 
%
%                   --> Optimized computation of Hermite polynomials in HERM_EVAL and HERM_EVAL_MULTIDIM
%
% -> 2017, Jan. 31  --> fixed a bug in REDUCE_SPARSE_GRID which might have affected problems with more than 10
%                       dimensions (i.e., some points that should have detected as equal might have considered
%                       different, causing some larger counting of sparse grid points and extra unnecessary
%                       work to evaluate functions)
%
%                   --> TENSOR_GRID returns also the set of points in each direction as a cell array (S.knots_per_dim) and the vector m
%                       (which is an input of TENSOR_GRID and denotes the number of points in each direction, m(n)==length(S.knots_per_dim{n}). 
%                       This information is carried over to sparse grids generated by SMOLYAK_GRID and
%                       SMOLYAK_GRID_MULTIIDX_SET, and then used to speed up INTERPOLATE_ON_SPARSE_GRID
%
%                   --> using LAGR_EVAL_FAST and further optimization of INTERPOLATE_ON_SPARSE_GRID for
%                       additional speedup
%
%                   --> added functions ISTENSOR and ISSMOLYAK to verify if a variable is a tensor / sparse
%                       grid
%
%                   --> added support for lexicographic operations: ISLEXICO compares two vectors and
%                       determines whether they are in lexicographic order; FIND_LEXICOGRAPHIC find a row in
%                       a matrix exploiting the fact that the matrix is lexicographically sorted.
%
%                   --> SMOLYAK_GRID and SMOLYAK_GRID_MULTIIDX_SET now accept as input a sparse grid, from
%                       which they will recover common tensor grids instead of recomputing them
%
%                   --> clarified help of REDUCE_SPARSE_GRID
%
%                   --> added error identifiers to support error handling (try/catch, rethrow). See ERRORS_AND_WARNINGS_LIST for a list 
%                       
% -> 2016, Jun. 20  --> fixing a bug in output computation in quadrature_on_sparse_grids.
% 
% -> 2015, Dec. 11  --> added function to ask number of open parallel matlab workers. Fixed a minor bug in close_parallel()
%
%
%
% RELEASE NOTES version 15.8 (Woodstock)
%
%
% -> 2015, Jun. 15  --> the new matlab toolbox syntax is now supported (see functions ACTIVATE_PARALLEL, CHECK_IF_PARALLEL_ON, CLOSE_PARALLEL)
%
%                   --> added MULTIIDX_GEN and MULTIIDX_BOX_SET to the tutorial
%
% -> 2015, Apr. 27  --> added function AS_REDUCED that converts a list of points / weights in a reduced grid
%                       structure
% 
%
%
% RELEASE NOTES version 14.12 (Fenice)
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
% RELEASE NOTES version 14.4 (Ritchie)
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
