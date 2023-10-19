function [O] = obsvt(A, C, T)
%OBSVT obsv for a set number of time steps forward
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 18, 2023

O = zeros(T*size(C,1), size(A,2));
Ai = eye(size(A));
for t=1:T
    O((t-1)*size(C,1)+1:t*size(C,1),:) = C*Ai;
    Ai = A * Ai;
end

end

