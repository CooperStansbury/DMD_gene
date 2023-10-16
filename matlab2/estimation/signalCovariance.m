function [R] = signalCovariance(signalsTrain)
%SIGNACOVARIANCE
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

% Remove the first time point from the data
signals = [];
for i=1:numel(signalsTrain)
    S = signalsTrain{i}; S(isnan(S)) = 0;
    signals = [signals; S(2:end,:)];
end

% least squares problem
R = cov(signals);

end