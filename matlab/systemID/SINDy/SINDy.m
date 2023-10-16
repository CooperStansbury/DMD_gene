function [eps] = SINDy(D,order,reps,cont)
%SINDY Fits a polynomial to time series data
%
% INPUTS
%   - D is n x t time series data (n) is features, (t) is time
%   - order: is polynomial order
%   - reps: number of replicates
%   - cont: if true is continuous, otherwise it is discrete
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 2, 2023

%% Argument handeling
if nargin == 3
    cont = false;
end

[n,T] = size(D);
t = T / reps;

%% Construct shifted data: FF and Xp

% If: the system is continuous
if cont
    F = firstDifferences(D, reps);
    FF = [];
    Xp = [];
    for i=1:reps
        FF = [FF F{i}];
        Xp = [Xp D(:,(i-1)*t + 1:i*t-1)];
    end
% Else: the system is discrete
else
    FF = [];
    Xp = [];
    for i=1:reps
        FF = [FF D(:,(i-1)*t + 2:i*t)];
        Xp = [Xp D(:,(i-1)*t + 1:i*t-1)];
    end
end

%% Construct theta
lib = makePolyLib(n, order);
theta = makeTheta(Xp, lib);

%% Solve least squares problem
% eps = theta' \ FF';

Xi = sparsifyDynamics(theta',FF',0.1,size(D,1));

end

