function isr = isreduced(S)

% ISREDUCED(S) returns 1 if S is a reduced sparse grid. A reduced sparse grid is a struct with fields 
% 'knots','m','weights','n'.


%----------------------------------------------------
% Sparse Grid Matlab Kit
% Copyright (c) 2009-2015 L. Tamellini, F. Nobile
% See LICENSE.txt for license
%----------------------------------------------------


if isstruct(S)
    isr=isempty(setxor(fieldnames(S),{'knots','m','weights','n'})) && length(S)==1;
else
    isr=0;
end


