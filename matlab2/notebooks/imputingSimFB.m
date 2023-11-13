%% Imputing Simulated Time Series Data with Forward Backward Dynamic Mode Decomposition
%
%   
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 13, 2023

close all; clc; clear;

%% Make Af and Ab

n = 2;
Af = [-0.001, -0.5; 0.5, -0.001];
Ab = [-0.001, 2; -2, -0.001];

T = 10000;
sfreq = 300;
delta = 0.01;
tt = delta * (1:T);

figure;
x0 = rand(2,1);
x0 = x0 / sum(x0);
tiledlayout(3,1);

% Simulate Af
X = zeros(n, T);
X(:,1) = x0;
for t=1:T-1
    X(:,t+1) = X(:,t) + delta * Af * X(:,t);
end
nexttile; plot(tt, X'); xlabel('Time'); ylabel('Amplitude'); title('Forward');

Abi = inv(Af);
X = zeros(n, T);
X(:,1) = x0;
for t=1:T-1
    X(:,t+1) = X(:,t) + delta * Abi * X(:,t);
end
nexttile; plot(tt, X'); xlabel('Time'); ylabel('Amplitude'); title('Backward');

X = zeros(n, T);
X(:,1) = x0;
for t=1:T-1
    c = t / T - 0.25;
    X(:,t+1) = X(:,t) + delta * ((1-c) * Af + c * Abi) * X(:,t);
end
nexttile; plot(tt, X'); xlabel('Time'); ylabel('Amplitude'); title('Mixed'); hold on;

X_sampled = X(:, 1:sfreq:end);
tt_sampled = tt(1:sfreq:end);

scatter(tt_sampled, X_sampled, 'k');



%%

theta = 30
Ab = [cos(theta) -sin(theta);
      sin(theta)  cos(theta)]

eig(Af)
eig(inv(Af))

eig(Ab)
eig(inv(Ab))

