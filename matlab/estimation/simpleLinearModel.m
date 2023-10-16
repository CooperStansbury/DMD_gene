function [lm] = simpleLinearModel(signalsTrain, velocityTrain)
%SIMPLELINEARMODEL COnstructs a simple linear model mapping PIP-FUCCI to
% velocity
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

% Remove the first time point from the data
signals = [];
vel = [];
for i=1:numel(signalsTrain)
    S = signalsTrain{i}; S(isnan(S)) = 0;
    V = velocityTrain{i};
    signals = [signals; S(2:end,:)];
    vel = [vel; V(2:end,:)];
end

% least squares problem
disp(size(signals));
lm = signals \ vel;

end