function isr = isreduced(S)

% ISREDUCED(S) returns 1 if S is a reduced sparse grid. A reduced sparse grid is a struct with fields 
% 'knots','m','weights','n'.

if isstruct(S)
    isr=isempty(setxor(fieldnames(S),{'knots','m','weights','n'})) && length(S)==1;
else
    isr=0;
end


