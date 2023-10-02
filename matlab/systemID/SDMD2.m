function [A] = SDMD2(DataMatrix, reps, S)
%SDMD Structured Dynamic Mode Decomposition
%
%   min_A |Xp - X * A| s.t. (A ~= 0) == S
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 14, 2023

%% Construct data matrices note, time snapshots are columns
t = size(DataMatrix,2) / replicates;
n = size(DataMatrix,1);
Xp = [];
Xf = [];
for r=1:replicates
    R = DataMatrix(:,t*(r-1)+1:t*r);
    Xp = [Xp R(:,1:end-1)];
    Xf = [Xf R(:,2:end)];
end

%% Least squares solution
A = Xp' \ Xf';

%% Sparsify Dynamics
% S is our sparsity strucutre and lambda zeros small parameters
for itr=1:maxItrs
    smallinds = (abs(A)<lambda);  % find small coefficients
    smallinds = smallinds + S;     % sparse structure
    A(smallinds > 0)=0;                % and threshold
    for ind = 1:n                   % n is state dimension
        biginds = ~smallinds(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        A(biginds,ind) = Xp(:,biginds)\Xf(:,ind); 
    end
end

end

