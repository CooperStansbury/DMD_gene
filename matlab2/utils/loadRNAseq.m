function [D, G] = loadRNAseq(ds, targetGenes)
%LOADRNASEQ
%
% ARGS:
%   - ds: 1 is 2015, 2 is 2017 data
%   - targetGenes: string array of genes of interest
%   - normalize:
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 13, 2023

if ds == 1
    [Dm,G,reps] = load2015(false); % Load data set
elseif ds == 2
    [Dm,G,reps] = loadMYOD(); % Load data set
end

if nargin == 1
    geneIdxs = 1:size(Dm,1);
else
    % subset the data to only include sensor and target genes
    gene2idx = containers.Map;                              % Map gene names to indices
    for i=1:numel(G); gene2idx(string(G{i})) = i; end
    geneIdxs = [];                                     % Identify cell cycle indices from the data
    for i=1:length(targetGenes)
        if gene2idx.isKey(string(targetGenes(i)))
            geneIdxs(end+1) = gene2idx(string(targetGenes(i)));
        end
    end
end

T = size(Dm,2)/reps;
D = zeros(numel(geneIdxs), size(Dm,2)/reps, reps);
for i=1:reps
    D(:,:,i) = Dm(geneIdxs,(i-1)*T+1:i*T);
end

G = G(geneIdxs);

end

