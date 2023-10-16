function [C] = getC(S, n)
%GETC Constructs observation matrix
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 18, 2023

m = numel(S);
C = sparse(m,n);
for i=1:m
    C(i,S(i)) = 1;
end

end

