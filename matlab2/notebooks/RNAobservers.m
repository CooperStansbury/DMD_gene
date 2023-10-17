%% RNA Observers
%
%   This file uses DMD + LQE/Kalman estimation to estimate gene expression
%   from a subset of genes
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 17, 2023

clear; close all; clc;
isPlot = 0;

ds = 2;
sensorGenes = ["PCNA","CDT1","GEM"];
targetGenes = [];
% Ovservable genes from CDT1, PCNA, and GEM according to PIP-FUCCI.m
targetGenes = ["PPP2R5B","CCND1","CCND2","CCND3","CDK4","CDK6","RB1","RBL1","RBL2","ABL1","HDAC1","HDAC2"];%,"E2F1","E2F2","E2F3","E2F4","E2F5","TFDP1","TFDP2","GSK3B","TGFB1","TGFB3","SMAD2","SMAD3","SMAD4"];

if ds == 1
    [D,G,reps] = load2015(false); % Load data set
elseif ds == 2
    [D,G,reps] = loadMYOD(); % Load data set
end

% subset the data to only include sensor and target genes
gene2idx = containers.Map;                              % Map gene names to indices
for i=1:numel(G); gene2idx(string(G{i})) = i; end
geneIdxs = [];                                     % Identify cell cycle indices from the data
for i=1:length(sensorGenes)
    if gene2idx.isKey(string(sensorGenes(i)))
        geneIdxs(end+1) = gene2idx(string(sensorGenes(i)));
    end
end
for i=1:length(targetGenes)
    if gene2idx.isKey(string(targetGenes(i)))
        geneIdxs(end+1) = gene2idx(string(targetGenes(i)));
    end
end
Dwhole = D;
D = D(geneIdxs,:);
G = G(geneIdxs);

% Construct system matrices (A,C)
out = shiftedDMD(D, reps, [], 1);
A = out.DMD.A_bar;
C = getC(1:3, length(geneIdxs)); C = full(C);

if isPlot
    E = eig(A); r = real(E); i = imag(E); theta = linspace(0, 2*pi, 1000); x = cos(theta); y = sin(theta);
    figure; plot(x, y); hold on; scatter(r,i); xlabel('Real ($\lambda$)','Interpreter', 'latex');
    ylabel('Imaginary ($\lambda$)','Interpreter', 'latex'); title('DMD Eigenvalues','Interpreter', 'latex');
end

%% Observability over the 3 replicates from 2018 MYOD
itrs = 1000;

B0 = zeros(size(A,1),1);    % Construct kalman filter
D0 = zeros(size(C,1),1);
sym = ss(A,B0,C,D0);
[kalmf,~,~] = kalman(sym,1e-4,1e-4);

nSensors = numel(sensorGenes);
T = size(D,2) / reps;
tspan = 1:T-1;
E = cell(3,1);
for r=1:reps
    xtrue = D(:,T*(r-1)+2:T*r);
    ytrue = xtrue(1:nSensors,:);
    ER = [];
    for i=1:itrs
        x0 = xtrue(:,1); % Set initial entry
        ind_1 = randi([1 size(Dwhole,1)],[numel(x0) - (nSensors), 1]);
        x0(nSensors + 1:end) = Dwhole(ind_1,2);
        % x0(nSensors + 1:end,:) = rand(numel(x0) - nSensors,1);
        xhat = lsim(kalmf,ytrue,tspan,x0); xhat = xhat';
        xhat = xhat(length(sensorGenes) + 1:end,:);
        ei = [];
        for t=1:size(xtrue,2)
            xh = xhat(:,t) / sum(xhat(:,t));
            xt = xtrue(:,t) / sum(xtrue(:,t));
            e = norm(xh - xt);   % we could use KL
            ei = [ei e];
        end
        ER = [ER; ei];
    end
    E{r} = ER;
end

% Random guessing
ER = [];
for i=1:itrs
    ei = [];
    for t=1:size(xtrue,2)
        % xh = xhat(:,t) / sum(xhat(:,t));
        xh = rand(); xh = xh / sum(xh);
        xt = xtrue(:,t) / sum(xtrue(:,t));
        e = norm(xh - xt);   % we could use KL
        ei = [ei e];
    end
    ER = [ER; ei];
end
E{r+1} = ER;

% Visualize the errors
figure;
tspan = 8*(1:size(ER,2));
for r=1:reps+1
    ER = E{r};
    m = mean(ER);
    err = std(ER) / sqrt(itrs);
    errorbar(tspan,m,err);
    if r == 1; hold on; end
end
title('LQE/Kalman Filter Estimator Error','Interpreter', 'latex');
legend(["Replicate 1", "Replicate 2", "Replicate 3","Random"]);
xlabel("Time (hours)",'Interpreter', 'latex')
ylabel("Norm between Estimated and True State",'Interpreter', 'latex');



%% state space model

B0 = zeros(size(A,1),1);    % Construct kalman filter
D0 = zeros(size(C,1),1);
sym = ss(A,B0,C,D0);
[kalmf,~,~] = kalman(sym,1e-4,1e-4);

% kalmf.OutputGroup
% kalmf.A

% Apply the kalman filter
tspan = 1:1:15;
x0 = D(:,2); x0(length(sensorGenes) + 1:end,:) = rand(numel(x0) - numel(sensorGenes),1);
y = lsim(kalmf,zeros(numel(tspan),3),tspan,x0); y = y';
xhat = y(length(sensorGenes) + 1:end,:);
xtrue = D(:,2:16);

% Error analysis
E = [];
for t=1:size(xtrue)
    e = norm(xhat(:,t) - xtrue(:,t));
    E = [E e];
end
figure; plot(E);


%%
[kalmf,L,P] = kalman(sym,0,0)
kalman(sym)
% Q = cov(D'); R = cov(D(1:length(sensorGenes),:)');