function [A,C,G,Q] = systemModel(ds)
%SYSTEMMODEL Linear system model for PIP-FUCCI and cell movement
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

sensorGenes = ["PCNA","CDT1","GEM"];
motileGenes = ["ROCK1", "ROCK2"];
% ts = 10;    % time scale: images taken once per 10 minutes. This is
% accounted for manually below when A is normalized


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

Q = cov(D');

end

