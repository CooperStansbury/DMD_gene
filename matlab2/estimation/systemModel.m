function [A,C,G,Q,D] = systemModel(ds, sensors, targetGenes, normalizerGene, timeshift)
%SYSTEMMODEL Linear system model for PIP-FUCCI and cell movement
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

sensorGenes = sensors; % ["PCNA","CDT1","GEM"];
% targetGenes = genes;   % ["ROCK1", "ROCK2"];
% ts = 10;    % time scale: images taken once per 10 minutes. This is
% accounted for manually below when A is normalized


if ds == 1
    [D,G,reps] = load2015(true); % Load data set
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
D = D(geneIdxs,:);
G = G(geneIdxs);

% normalize the data according to the normalizerGene
%{
for i=1:length(G); if strcmp(G{i}, normalizerGene); idxN = i; end; end
t = size(D,2)/reps;
for i=1:reps
    repI = D(:,(i-1)*t+1:i*t);
    [~, locOfMaxExpressionOfNormalizerGene] = max(repI(idxN,:));
end
%}

% Construct system matrices (A,C)
out = shiftedDMD(D, reps, [], 1);
A = out.DMD.A_bar;
% Mapping between time scales
%       A = exp(A' * 8 hours) --> log(A) / (8*60 min) = A'
%       B = exp(A' * 10 min)  --> B = exp(log(A) * 10 min / (8 * 60 min))
if timeshift ~= -1
    A = exp(log(A) * timeshift / 480);
end
C = getC(1:3, length(geneIdxs)); C = full(C);

Q = cov(D');

end

