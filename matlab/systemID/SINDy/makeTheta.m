function theta = makeTheta(D, lib)
%MAKETHETA
%
%   INPUTS:
%       - D: time series data matrix is features (n) by time (t)
%       - lib: library of nonlinear terms (p) to be placed in theta
%   OUTPUTS:
%       - theta: nonlinear functions by time is a nonlinear term (p) by
%       time (t) matrix
%
%{
    D = rand(3,10);
    F = firstDifferences(D,1); F = F{1};
    lib = { @(D)D(1,:),
            @(D)D(2,:),
            @(D)D(3,:),
            @(D)D(1,:).*D(2,:),
            @(D)D(1,:).*D(3,:),
            @(D)D(2,:).*D(3,:)  };
%}
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 2, 2023

[~,t] = size(D);
p = length(lib);
theta = zeros(p,t);

for i=1:p
    f = lib{i};
    theta(i,:) = f(D);
end

end