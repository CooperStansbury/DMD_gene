%% Observers
%
%   In this file we apply dynamic observers on PIP-FUCCI sensors to make
%   predictions about hidden or unobserved cell states.
%
%   System Output/Sensors: The PIP-FUCCI assay generates three channel
%   images where each channel may be utilized to estimate the expression of
%   the CDT1, PCNA, and GEM genes. We are not able to directly determing
%   the gene expression similar to RNAseq but rather know the luminance of
%   different color channel correlates with gene expression
%
%   Validation:
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

%% Predicting Velocity
clear; close all; clc;

% load system model
ds = 1;
[A,C,G,Q] = systemModel(ds);

% Load PF signals
minSamples = 10;
cellTracks = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
cellTracks = filterCellTracks(cellTracks,minSamples);
velocity = getVelocity(cellTracks);
signals = getPIPFUCCI(cellTracks);
Q = signalCovariance(cellTracks);

% Partition data for testing
kfold = 5;
nsamples = numel(signals);
cv = cvpartition(nsamples, 'kfold', kfold);

%% Perform cross fold validation
nskip = 5;
for fold=1:kfold
    % Split data for a fold
    idxs = training(cv,fold);
    train = find(idxs == 1);
    test = find(idxs == 0);
    signalsTrain = signals(train);
    signalsTest = signals(test);
    velocityTrain = velocity(train);
    velocityTest = velocity(test);
    R = signalCovariance(signalsTrain);

    % Construct linear models on training data
    slm = simpleLinearModel(signalsTrain, velocityTrain);
    klm = kalmanLinearModel(signalsTrain, velocityTrain, A, C, nskip, Q, R);

    % Apply linear models on testing data

end


unique([train; test])
testing(cv,1)


%% Sandbox Below

%% Build Linear System Model of PIP-FUCCI Sensors with Motility Genes
sensorGenes = ["PCNA","CDT1","GEM"];
motileGenes = ["ROCK1", "ROCK2"];
ds = 1;
ts = 10;    % time scale: images taken once per 10 minutes


if ds == 1
    [D,G,reps] = load2015(true); % Load data set
elseif ds == 2
    [D,G,reps] = loadMYOD(); % Load data set
end

gene2idx = containers.Map;                              % Map gene names to indices
for i=1:numel(G); gene2idx(string(G{i})) = i; end
geneIdxs = [];                                     % Identify cell cycle indices from the data
for i=1:length(sensorGenes)
    if gene2idx.isKey(string(sensorGenes(i)))
        geneIdxs(end+1) = gene2idx(string(sensorGenes(i)));
    end
end
for i=1:length(motileGenes)
    if gene2idx.isKey(string(motileGenes(i)))
        geneIdxs(end+1) = gene2idx(string(motileGenes(i)));
    end
end

D = D(geneIdxs,:);
G = G(geneIdxs);

% Construct system matrices (A,C)
out = shiftedDMD(D, reps, [], 0.9);
A = out.DMD.A_bar;
% Mapping between time scales
%       A = exp(A' * 8 hours) --> log(A) / (8*60 min) = A'
%       B = exp(A' * 10 min)  --> B = exp(log(A) * 10 min / (8 * 60 min))
A = exp(log(A) * 10 / 480);
C = getC(1:3, length(geneIdxs)); C = full(C);

%% Extract System Output (Y) and Function/State (speed)
MV = readtable('data/C1.tracks.full.csv'); % Read in cell tracks

% Construct system output
systemOutput = [MV.ID MV.c0_intensity_mean MV.c1_intensity_mean MV.c2_intensity_mean];
id = mode(systemOutput(:,1)); % select the cell with the most samples
systemOutput(systemOutput(:,1) ~= id,:) = [];
Y = systemOutput(:,2:4);

% Construct cell velocity 
cellMeasurements = MV;
cellMeasurements(cellMeasurements.ID ~= id,:) = [];
velocity = zeros(height(cellMeasurements),1);
for t=2:numel(velocity)
    % velocity = distance / time
    velocity(t) = sqrt((cellMeasurements.x(t) - cellMeasurements.x(t-1))^2 + ...
        (cellMeasurements.y(t) - cellMeasurements.y(t-1))^2) / ts;
end

%% Perform Kalman Filtering

n = size(A,1);
Y = Y';

% Number of time steps
T = size(Y, 2);

% Initialize Kalman filter parameters
% x_hat_0 = rand(n,1); % Replace with your initial estimate
x_hat_0 = [Y(:,1); mean(Y(:,1)); mean(Y(:,1))];
P_0     = zeros(n,n);   % Replace with your initial covariance matrix
Q       = cov(D');   % Replace with your process noise covariance
R       = cov(Y');    % Replace with your measurement noise covariance

% Initialize variables
x_hat = x_hat_0;
P = P_0;
% Y % = % Y';

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

figure; plot(real(Xhat))

%% Linear Models to relate sensor output with velocity
lm = Y' \ velocity;     % lm maps system output directly to velocity
Xhat(:,1:3) = Y';
km = Xhat \ velocity;   % km maps filtered system output to velocity

%% Plot error over time for each map
lerror = lm' * Y;
kerror = km' * Xhat';

figure;
scatter(1:T, lerror); hold on;
scatter(1:T, kerror);

mn = 12; figure;
mmk = movmean(kerror, mn); mmk = real(mmk);
mml = movmean(lerror, mn);
plot(mml); hold on;
plot(mmk)
legend(["Linear Model","Kalman Filtering"]);


%% How to deal with inevenly spaced RNAseq and PIP-FUCCI Samples
%       A = exp(A' * 8 hours) --> log(A) / (8*60 min) = A'
%       B = exp(A' * 10 min)  --> B = exp(log(A) * 10 min / (8 * 60 min))
