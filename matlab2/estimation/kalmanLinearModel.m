function [lm] = kalmanLinearModel(signalsTrain, velocityTrain, A, C, nfirst, Q, R)
%KALMANLINEARMODEL Constructs a linear model mapping the Kalman filtering 
% state estimation from PIP-FUCCI signals to velocity
%
%   INPUTS:
%       - signalsTrain: training data
%       - velocityTrain: training labels
%       - A: linear model
%       - C: system inputs
%       - nfirst: number of time steps to ignore while the filter converges
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

n = size(A,1);

% set kalman filter parameters
P_0     = eye(n,n);   % Replace with your initial covariance matrix

% kalman filter each signal
signals = [];
vel = [];
for i=1:numel(signalsTrain)
    % Extract individual signal
    S = signalsTrain{i}; S(isnan(S)) = 0;
    Y = S; % rename variable so it fits the Kalman filtering code
    V = velocityTrain{i};

    T = size(Y, 1); % Number of time steps
    Y = Y';
    % Initialize Kalman filter parameters
    x_hat_0 = [Y(:,1); mean(Y(:,1)); mean(Y(:,1))];
    x_hat = x_hat_0;
    P = P_0;

    % Kalman filter loop
    Xhat = zeros(T, n);
    for t = 1:T
        % Prediction step
        x_hat_minus = A * x_hat;            % Predicted state estimate
        P_minus = A * P * A' + Q;           % Predicted error covariance
        
        % Update step (using the measurement Y(:, t))
        K = P_minus * C' / (C * P_minus * C' + R); % Kalman gain
        x_hat = x_hat_minus + K * (Y(:, t) - C * x_hat_minus); % Updated state estimate
        P = (eye(size(A)) - K * C) * P_minus; % Updated error covariance
        Xhat(t,:) = x_hat;
    end

    % save data
    signals = [signals; real(Xhat(nfirst:end,:))];
    vel = [vel; V(nfirst:end,:)];
end

% least squares problem
disp(size(signals));
lm = signals \ vel;

end