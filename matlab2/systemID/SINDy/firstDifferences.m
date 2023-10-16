function [F] = firstDifferences(D, reps)
%FIRSTDIFFERENCES Calculates the first difference or discrete derivative
%from time series samples.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 2, 2023

t = size(D,2) / reps;
F = cell(reps,1);
for r=1:reps
    Dr = D(:,(r-1)*t+1:r*t);
    Fr = zeros(size(Dr,1), size(Dr,2)-1);
    for i=1:t-1
        Fr(:,i) = Dr(:,i+1) - Dr(:,i);
    end
    F{r} = Fr;
end

end

