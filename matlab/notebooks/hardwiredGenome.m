%% Hardwired Genome
%
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 24, 2023

clear
HWG = load_HWG();


%% DMD of time series data

[D,G] = loadMYOD();             % Load MYOD data set
% A = shiftedDMD(D,3,[],0.9);
A = EDMD(D);
help dmdsp

%% Align DMD with Hardwired Genome

% Step 1: Get the indices of the k largest entries in A
k = 100;
[~, linIdx] = maxk(A(:), k);
[x,y] = ind2sub(size(A), linIdx);
I = [x, y];

% Step 2: Map indices to gene names from time series data
genePairs = cell(size(I,1), 2);
for i=1:size(genePairs,1)
    genePairs{i,1} = G(I(i,1));
    genePairs{i,2} = G(I(i,2));
end

% Step 3: look up gene pairs in the hardwired genome
ensg = HWG.geneIndexTable.("Stable ID");
ensg2idx = containers.Map;
for i=1:numel(ensg)
    if ~isempty(ensg{i})
        ensg2idx(ensg{i}) = i;
    end
end

gene2idx = containers.Map(HWG.geneIndexTable.("Gene Name"), 1:height(HWG.geneIndexTable));
gene2idx = containers.Map;
for i=1:numel(gn)
    if ~isempty(gn{i})
        gene2idx(gn{i}) = i;
    end
end

gene2idx.keys
gn = HWG.geneIndexTable.("Gene Name");
gene2idx = containers.Map(gn{:}, 1:height(HWG.geneIndexTable));

%% trash code
tic; disp(maxk(A(:), 10)); toc
k=1000;

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
