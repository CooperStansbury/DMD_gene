function [DC,GC,replicates] = loadTimeSeriesCC(ds)
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 3, 2023

cellCycleGenes = readtable('data/kegg_hsa04110.csv');   % Load cell cycle genes
if ds == 1
    [D,G,replicates] = load2015(true); % Load data set
else
    [D,G,replicates] = loadMYOD(); % Load data set
end
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

end