% ============================================================
%  Generates a tensor grid with m=[m1,m2,...,m_N] Gauss points 
%  in each direction and computes the corresponding weights
%
%  The information is stored in a three field structure S:
%    S.knots: vector containing the tensor grid knots
%    S.weights: vector containing the corresponding weights
%    S.size: size of the tensor grid = prod(m)
%
%  usage:
%   [S] = tensor_grid(N,m,knots)
%   input
%     N: dimension
%     m: vector containing the number of knots in each direction
%     knots: function defining the 1D gauss knots
%            header: [x,w]=knots(m)  (generates m knots and weights)
%   output
%     S: structure containing the information on the tensor grid (see above)
% ============================================================

function [S] = tensor_grid(N,m,knots)

[S.knots,S.weights] = knots(m(1));   
for i=2:N
 [Xi,Wi] = knots(m(i));
 S.knots = combvec(S.knots,Xi);         
 S.weights = combvec(S.weights,Wi);
 S.weights = prod(S.weights);
end

S.size = prod(m);
