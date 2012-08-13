% DIFFERENCES W.R.T. THE PREVIOUS RELEASE : 
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
%   admissiilbility condition
%
% -> /tools/polynomials now includes code to generate Legendre and Hermite polynomials. 
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
%     to rescale knots and weights after having generated the grid. See the help for more detailss
%         
% -> however, for more complex cases (e.g. each rand var lives on its own interval), one still needs to
%    rescale afterwards. To this end, the function  get_interval_map.m      
%    will help generating the function that rescales the knots and the quadrature weights
% 
% -> knots_trap and midpoint have been removed
% 
% -> convert_to_modal.m function now allows to compute the gPCE of a function starting from its sparse grid approximation.
%     The alogrithm implemented recasts the sum of lagrange polynomials into a sum of orthornormal polynomials and has 
%     been found to be more efficient than computing the gPCE coefficients by sparse grid quadrature.
%     Works for both Legendre and Hermite gPCE 