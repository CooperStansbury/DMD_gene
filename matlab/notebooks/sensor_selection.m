%% Sensor Selection
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 17, 2023

clear; close all; clc

%% 2018 MYOD Data Set
T = readtable("time_series_rna_fold_changes.csv");  % Read .csv file
fid = fopen("time_series_rna_fold_changes.csv",'r');
genes = fgetl(fid); genes = split(genes, ',');
genes = genes(2:end);
D1 = zeros(15,size(T,2)-1);                         % Create matrices to store data
D2 = zeros(15,size(T,2)-1);
D3 = zeros(15,size(T,2)-1);
for t=1:15                                          % Read each timepoint
    D1(t,:) = T{3*(t-1) + 1, 2:size(T,2)};
    D2(t,:) = T{3*(t-1) + 2, 2:size(T,2)};
    D3(t,:) = T{3*(t-1) + 3, 2:size(T,2)};
end

D1 = D1'; D2 = D2'; D3 = D3';
D = [D1 D2 D3];                                     % Organize the data
clearvars -except D genes                           % Remove all other variables

% Perform sensor selection
[S, GV] = hasnain2023(D);
genes(S)

%% 2023 Hasnain Data set

filename = "hasnain2023_data_fc_gene_list.csv";
genes = readtable(filename,'Delimiter','\t');

filename = "hasnain2023_data_fc.csv";
T = readtable(filename);
D = T{2:end,2:end};
clearvars -except D genes                           % Remove all other variables

% Perform sensor selection
[S, GV] = hasnain2023(D);
genes.GeneName(S)

%% Kalman Filtering
%
%   In this section I construct Kalman filters to estimate the full state
%   vector of gene expression. We evaluate the performance of each gene set
%   according to the convergence of each estimator.

clear
[D, genes] = loadMYOD(); replicates = 3; t = 15;
[D, genes] = loadHasnain2023(); replicates = 2; t = 9;
n = size(D,1);

output = shiftedDMD(D,replicates,[],0.9);

A_bar = output.Xp*output.DMD.VX*pinv(output.DMD.Sig)*output.DMD.UX';

% Perform sensor selection
[S, GV] = hasnain2023(D);
C = getC(S, n); C = full(C);
m = size(C,1);

Y = C * D; Y = Y(:,1:t);

% Vd = 1e-3 * ones(n,n); % process covariance
% Vn = 1e-3 * ones(m,m); % sensor covariance
% Vd = corr(D',D');
% Vn = corr(Y',Y');
Vd = 1e-1 * eye(n,n);
Vn = 1e-1 * eye(m,m);
B = zeros(n,1);
[L, ~, ~] = lqe(A_bar, Vd, C, Vd, Vn);

% Kf = lqr(A',C',Vd,Vn)';
sysKF = ss(A - L*C, [B L], eye(n,n), 0 * [B L]);





% Construct the system and filter
Ts = -1;    % Discrete system

sys = ss(A_bar,zeros(n,1),C,zeros(m,1),Ts,'OutputName','y');    % Construct system

Q = 0;    % process noise // We will need to determine how to estimate these values
R = 0;      % sensor noise

[kalmf,L,~,Mx,Z] = kalman(sys,Q,R); % Construct Kalman Filter

Xhat = lsim(kalmf, Ydata, t=1:9)
