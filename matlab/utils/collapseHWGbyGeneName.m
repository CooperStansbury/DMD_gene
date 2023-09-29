function [A, G] = collapseHWGbyGeneName(HWG)
%COLLAPSEHWGBYGENENAME Collapses HWG A matrix so that each row corresponds
% to a unique gene name
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 26, 2023

gn = HWG.geneIndexTable.("Gene Name");
gene2idx = containers.Map;
for i=1:numel(gn)
    if ~isempty(gn{i})
        gene2idx(gn{i}) = i;
    end
end

np = size(gene2idx.keys);
Ap = zeros(np, np);



end

