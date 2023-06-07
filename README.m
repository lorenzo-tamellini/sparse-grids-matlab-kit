%----------------------------------------------------------------------------------
% Sparse Grids Matlab Kit
% Copyright (c) 2009-2023 L. Tamellini, F. Nobile, C. Piazzola
% See LICENSE.txt for license
%----------------------------------------------------------------------------------



%---------------------------
% Contributors
%---------------------------
%
% 1) Lorenzo Tamellini
% 2) Diane Guignard
% 3) Fabio Nobile
% 4) Chiara Piazzola
% 5) Giovanni Porta
% 6) Bjorn Sprungk
% 7) Francesco Tesei


%---------------------------
% 1) QUICK START / RESOURCES
%
% Take a look at:
%
% -- RELEASE_NOTE.m for a detailed list of diffrences with the previous versions
%
% -- SPARSE_GRIDS_TUTORIAL.m in the folder docs-examples for a gradual introduction to the code
%
% -- sparse grids matlab kit USER MANUAL, available at https://sites.google.com/view/sparse-grids-kit/home



%---------------------------
% 2) HOW TO INSTALL THE TOOLKIT
% 
% add to path this folder (with subfolders). Alternatively, run 
%
% addpath(genpath(pwd))
%
% at the beginning of the Matlab session. Consider launching Matlab as  
%
% matlab -r "addpath(genpath(pwd))"
%
% under Unix, it will run add_to_path during the start-up




%---------------------------
% 3) HOW TO USE THE TOOLKIT
%
% To create a sparse grids, use any of these functions:
%
% -> create_sparse_grid
% -> create_sparse_grid_quick_preset
% -> create_sparse_grid_multiidx_set
% -> create_sparse_grid_add_multiidx
% -> adapt_sparse_grid

% A number of functionalities to operate on sparse grids are provided: 
%
% -> compute_sobol_indices_from_sparse_grid
% -> convert_to_modal
% -> derive_sparse_grid
% -> evaluate_on_sparse_grid
% -> export_sparse_grid_to_file
% -> hessian_sparse_grid.m
% -> interpolate_on_sparse_grid
% -> plot_sparse_grid
% -> plot3_sparse_grid
% -> plot_sparse_grid_interpolant
% -> quadrature_on_sparse_grid
% -> reduce_sparse_grid
% 
% type HELP FUNCTIONAME to access help for each function.
% Documentation and examples can be found in docs-examples folder. 




%---------------------------
% 4) REPORT BUGS to tamellini@imati.cnr.it
% 
% also, let us know your email address if you want us to warn you whenever we release a new version of the toolbox



%---------------------------
% 5) PLEASE CITE US

% Please cite our toolbox by mentioning the webpage containing the package (https://sites.google.com/view/sparse-grids-kit/home) 
% and adding the following reference to your work (check website for the most updated citation data):
%
% @article{piazzola.tamellini:SGK,
%  author = {Piazzola, C. and Tamellini, L.},
%  title  = {{The Sparse Grids Matlab kit - a Matlab implementation of sparse grids for high-dimensional function approximation and uncertainty quantification}},
%  journal= {},
%  year   = {2023},
%  volume = {},
%  number = {},
%  pages = {},
%}




