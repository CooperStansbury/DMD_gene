%% PIP-FUCCI
%
%   We investigate the observability of the cell cycle according to DMD
%   applied to the key genes in the cell cycle.
%
%   PIP FUCCI NOTES:
%       Genes of Interest: CDT1, PCNA, GEM
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 29, 2023

%% Load data
clear; close all; clc;
cellCycleGenes = readtable('data/kegg_hsa04110.csv');   % Load cell cycle genes
[D,G] = loadMYOD(); replicates = 3;                     % Load data set
gene2idx = containers.Map;                              % Map gene names to indices
for i=1:numel(G); gene2idx(string(G{i})) = i; end
cellCycleIdxs = [];                                     % Identify cell cycle indices from the data
for i=1:height(cellCycleGenes)
    if gene2idx.isKey(string(cellCycleGenes.GeneName{i}))
        cellCycleIdxs(end+1) = gene2idx(string(cellCycleGenes.GeneName{i}));
    end
end

GC = G(cellCycleIdxs);      % Extract cell cycle gene indices
DC = D(cellCycleIdxs,:);    % Extract cell cycle RNAseq

clearvars -except GC DC replicates

%% Is the system observable from GEM, PCNA, and CDT1?
sensors = ["GEM", "PCNA", "CDT1"];
sensorIdxs = [];
for s=1:length(sensors)
    for i=1:length(GC)
        if strcmp(sensors(s), string(GC{i}))
            sensorIdxs(end+1) = i;
            break;
        end
    end
end
C = getC(sensorIdxs, size(DC,1)); C = full(C);
dmd = shiftedDMD(DC,3,[],0.9);
A = dmd.Xp*dmd.DMD.VX*pinv(dmd.DMD.Sig)*dmd.DMD.UX';
O = obsv(A,C);
disp(rank(O));
[U,sigma,V] = svd(O);

figure; scatter(1:min(size(sigma)), diag(sigma),'.'); set(gca,'YScale','log');
title('Observability Matrix Scree Plot'); ylabel('Singular Value'); xlabel('Component');
figure; imagesc(A); title('DMD Linear Operator');
figure; plot(1:size(DC,2),DC(sensorIdxs,:));
title("PIP-FUCCI Sensors"); xlabel("Time"); ylabel("Expression"); legend(sensors);

%% Perform Sensor Selection
[s, GV] = hasnain2023(DC, replicates);
GC(s(1:5))

%% Perform DMD
dmd = shiftedDMD(DC,3,[],0.9);



