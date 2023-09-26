%% Hardwired Genome
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 24, 2023

clear
HWG = load_HWG();

[D,G] = loadMYOD();
A = shiftedDMD(D,3,[],0.9);
A = EDMD(D);
k=1000;
% Step 1: Get the indices of the k largest entries in A
[sortedA, sortedIndices] = sort(A(:), 'descend');
kLargestIndices = sortedIndices(1:k);

ItoHWGixdMap = zeros(numel(G),1);
HWGG = HWG.geneIndexTable.("Gene Name");
for i=1:length(G)
    geneNameInA = G{i};
    rowIdxInHWGTable = find(strcmp(geneNameInA, HWGG));
    if isempty(rowIdxInHWGTable)
        rowIdxInHWGTable = -1;
    end
    ItoHWGixdMap(i) = rowIdxInHWGTable;
    disp(i);
end

% Step 2: Map the indices of A to the indices of HWG.A using the two lists of indices
rowIndicesInHWGA = cellfun(@(i) G{i}(1), num2cell(kLargestIndices));
colIndicesInHWGA = cellfun(@(i) G{i}(2), num2cell(kLargestIndices));


% [D,G] = loadHasnain2023();
