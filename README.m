%----------------------------------------------------------------------------------
% Sparse Grids Matlab Kit
% Copyright (c) 2009-2017 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------------------------------------


% Please cite our toolbox by mentioning the webpage containing the package (http://csqi.epfl.ch) 
% and adding the following reference to your work:
%
% @InCollection{back.nobile.eal:comparison,
%  author    = {B\"ack, J. and Nobile, F. and Tamellini, L. and Tempone, R.},
%  title        = {Stochastic spectral {G}alerkin and collocation methods for {PDE}s with random coefficients: a numerical comparison},
%  booktitle    = {Spectral and High Order Methods for Partial Differential Equations},
%  pages        = {43--62},
%  publisher    = {Springer},
%  year        = 2011,
%  volume    = 76,
%  series    = {Lecture Notes in Computational Science and Engineering}, 
%  editor    = {Hesthaven, J.S. and Ronquist, E.M.},
%  note        = {Selected papers from the ICOSAHOM '09 conference, June 22-26, Trondheim, Norway}
%}


%---------------------------
% Take a look at RELEASE_NOTE.m for a detailed list of diffrences with the previous version



%---------------------------
% 1) HOW TO INSTALL THE TOOLKIT
% 
% add to path this folder (with subfolders). Alternatively, run 
%
% addpath(genpath(pwd))
%
% at the beginning of the Matlab session. Consider launching Matlab as  
%
% matlab -r addpath(genpath(pwd))
%
% under Unix, it will run add_to_path during the start-up




%---------------------------
% 2) HOW TO USE THE TOOLKIT
%
% The functionalty provided are: 
%
% -> smolyak_grid.m
% -> smolyak_grid_multiidx_set.m
% -> reduce_sparse_grid.m
% -> evaluate_on_sparse_grid.m
% -> interpolate_on_sparse_grid.m
% -> quadrature_on_sparse_grid.m
% -> convert_to_modal.m
%
% documentation and examples can be found in docs-examples folder. 
%
% In particular, docs-examples/html contains the published version of sparse_grid_tutorial.m (file/publish).



