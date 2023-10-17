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
%   TODO:
%       1. modify systemModel to accept a time scale argument for how the
%       DMD and PF matrices work together
%
%   PIP-FUCCI: we just used fucci reporting system
%   - red = G1
%   - both = S
%   - green = S,G2,M
%
%   IMAGES 2 GENE EXPRESSION:
%       Problems:
%           1. we don't know the cell state of cells at the start of the
%           imaging, so fold change is bad because each cell sample will
%           not be aligned with one another
%           2. 
%
%   cellTracks
%   - c0: red (G1)
%   - c1: nucleus
%   - c2: green
%   - prob: probability of segmentation (0.4 was set previously by CS)
%   - label: just an id, not needed
%   - state != state in cell cycle
%   - dummy indicates if a cell is imputed in a time point or actually
%   detected in the image. If dummy is true, then the rest is NAN
%   - mean: dominated by poor segmentation or outliers
%   - foldChange: compares intensity of a cell to the intensity of the rest
%   of the image, use the processed ones
%   - pval: significance of the fold change based on number of samples
%   (pixels) in the image
%   - estimate of image/pixel intensity is performed as:
%       1. time series filtering
%       2. static estimate
%       3. smoothing of static estimates
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

%% Predicting Velocity
clear; close all; clc;

ds = 2;
% normalizerGene = "GEM"; % "PNCA" or "CDT1"
sensors = ["PCNA","CDT1","GEM"];
% targetGenes = ["ROCK1", "ROCK2"];
% Ovservable genes from CDT1, PCNA, and GEM according to PIP-FUCCI.m
targetGenes = ["PPP2R5B","CCND1","CCND2","CCND3","CDK4","CDK6","RB1","RBL1","RBL2","ABL1","HDAC1","HDAC2","E2F1","E2F2","E2F3","E2F4","E2F5","TFDP1","TFDP2","GSK3B","TGFB1","TGFB3","SMAD2","SMAD3","SMAD4"];

% load system model
[A,C,G,Q] = systemModel(ds, sensors, targetGenes); % 4th argument is ignored because Q covariance is calculated below from the cell tracks, not from

% Load PF signals
minSamples = 10;
cellTracks = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
cellTracks = filterCellTracks(cellTracks,minSamples, true);
velocity = getVelocity(cellTracks);
signals = getPIPFUCCI(cellTracks);
Q = signalCovariance(signals);

% Partition data for testing
kfold = 5;
nsamples = numel(signals);
cv = cvpartition(nsamples, 'kfold', kfold);

%% Perform cross fold validation
nskip = 5;
for fold=1:kfold
    % Split data for a fold
    idxs          = training(cv,fold);
    train         = find(idxs == 1);
    test          = find(idxs == 0);
    signalsTrain  = signals(train);
    signalsTest   = signals(test);
    velocityTrain = velocity(train);
    velocityTest  = velocity(test);
    R = signalCovariance(signalsTrain);

    % Construct linear models on training data
    slm = simpleLinearModel(signalsTrain, velocityTrain);
    klm = kalmanLinearModel(signalsTrain, velocityTrain, A, C, nskip, Q, R);

    % Apply linear models to predict velocity on testing data
    linearPredictions = predictLinear(slm, signalsTest);
    kalmanPredictions = predictKalman(klm, signalsTest, A, C, Q, R);

    % Evaluate predictions
    E = [];
    for j=1:numel(linearPredictions)
        l = linearPredictions{j}; l(isnan(l)) = 0;
        k = kalmanPredictions{j}; k(isnan(k)) = 0;
        v = velocityTest{j};

        l = l(end-3:end);
        k = k(end-3:end);
        v = v(end-3:end);

        le = sum(abs(v-l'));
        ke = sum(abs(v-k'));
        E = [E; le ke];
        % disp(string(le) + " v. " + string(ke));
        % figure;
        % plot(l); hold on;
        % plot(k);
        % plot(v);
        % legend(["Linear","Kalman","True"]);
        % title(string(le) + " v. " + string(ke));
    end
    mean(E)
    figure; plot(E)

end


unique([train; test])
testing(cv,1)


%% Figures

%% Visualize the convergence of the Kalman Filter


%% Eigenvalues of A

E = eig(A);
r = real(E);
i = imag(E);

% Define the angle values
theta = linspace(0, 2*pi, 1000);  % Create 1000 equally spaced points around the unit circle

% Compute the x and y coordinates of the unit circle
x = cos(theta);
y = sin(theta);

% Plot the unit circle
figure; plot(x, y); hold on;
scatter(r,i);
xlabel('Real ($\lambda$)','Interpreter', 'latex');
ylabel('Imaginary ($\lambda$)','Interpreter', 'latex');
title('DMD Eigenvalues','Interpreter', 'latex');

%%
% Luminiscinece and gene activity distribution
clear; close all; clc;

sensors = ["PCNA","CDT1","GEM"];
genes = ["ROCK1", "ROCK2"];

% load system model
ds = 1;
[~,~,G,~,D] = systemModel(ds, sensors, genes); % 4th argument is ignored because Q covariance is calculated below from the cell tracks, not from

D = D(1:3,:);
figure; scatter3(D(1,:),D(2,:),D(3,:),'.');

cellTracks = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
cellTracks = filterCellTracks(cellTracks, 10);
signals = getPIPFUCCI(cellTracks,unique(cellTracks.ID),'log','PCNA');
for i=1:numel(signals); Y = [Y; signals{i}]; end
figure; scatter3(Y(:,1),Y(:,2),Y(:,3),'.');


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
