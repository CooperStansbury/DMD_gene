function pred = predictKalman(klm, signalsTest, A, C, Q, R)
%PREDICTKALMAN Predict velocity with Kalman linear model
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

n = size(A,1);

filteredSignals = cell(numel(signalsTest), 1);

% set kalman filter parameters
P_0     = eye(n,n);   % Replace with your initial covariance matrix
for i=1:numel(signalsTest)
    S = signalsTest{i}; S(isnan(S)) = 0;
    Y = S; % rename variable so it fits the Kalman filtering code
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
    filteredSignals{i} = real(Xhat);
end

pred = predictLinear(klm, filteredSignals);

end
