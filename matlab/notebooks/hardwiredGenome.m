%% Hardwired Genome
%
%   We investigate the relationship between the full A operator identified
%   by DMD and the sparsity of the Hardwired genome.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 24, 2023

%% Build zero threshold hardwired genome
%
% Should run in about 30-40 minutes
clear all
close all
clc

indexTable = buildIndexTable();

% Preamble
thresh = 0; % Include all possible data

% Get adjacency lists
A_list_HuRI = list_HuRI();                              % get HuRI adjacency list
disp("    HuRI Data Accessed");
[A_list_STRING, ~] = list_STRING(thresh);               % get STRING adjacency list at specific threshold
disp("    STRING Data Accessed");
A_list = list_combine({A_list_HuRI, A_list_STRING});    % Combine the lists
% tic; A_HWG = list2mat(A_list, indexTable); disp(toc);   % build the matrix

unique_ids = indexTable.("Stable ID");
n = length(unique_ids);
gene_idxs = containers.Map;
for idx=1:length(unique_ids)
    gene_idxs(unique_ids(idx)) = idx;
end

tic;
E = zeros(height(A_list), 2);
for PPI=1:height(A_list)    % takes about 35 minutes to run
    if mod(PPI,1e6) == 0; disp(PPI); disp(toc); end
    if isKey(gene_idxs, A_list.("Protein 1")(PPI)) && isKey(gene_idxs, A_list.("Protein 2")(PPI))
        E(PPI,1) = gene_idxs(A_list.("Protein 1")(PPI));
        E(PPI,2) = gene_idxs(A_list.("Protein 2")(PPI));
        % A(gene_idxs(A_list.("Protein 1")(PPI)), gene_idxs(A_list.("Protein 2")(PPI))) = 1;
    end
end
disp(toc);
EE = E;
EE(all(EE == 0, 2), :) = [];
A = sparse(EE(:,1), EE(:,2), ones(size(EE,1), 1));

interacting_genes = find(sum(A) > 0);
A_HWG = A(interacting_genes, interacting_genes);
A_IndexTable = indexTable(interacting_genes, :);

% Construct the factor matrices and package as the Hardwired Genome
[B_HWG, C_HWG, C_IndexTable] = HWGMatrixDecomp(A_HWG, A_IndexTable);

% Construct and save HWG object
HWG = struct;
HWG.thresh = thresh;
HWG.A = A_HWG;
HWG.B = B_HWG;
HWG.C = C_HWG;
HWG.geneIndexTable = A_IndexTable;
HWG.TFIndexTable = C_IndexTable;

% open the correct directory to save this to
% save('HWG.mat', "HWG")
% nnz(HWG0.A) * 100 / numel(HWG0.A)
% save

%%

tic
A_list_HuRI = list_HuRI();
disp(toc); % ~4.866 seconds on JP's laptop
tic;
A_list_STRING = list_STRING(0);
disp(toc); % 1.951106600700000e+03 seconds on JP's laptop
tic;
A_lists = {A_list_HuRI, A_list_STRING};
A_list = list_combine(A_lists);
disp(toc); % 1.58 seconds
tic;
[A_HWG_0, HWG_0_info] = build_A(A_list);
disp(toc);

indexTable = buildIndexTable();

%% Sep. 27 2023

clear
% Select k "most important" gene-gene interactions identified by DMD
[D,G] = loadMYOD();
DMD = shiftedDMD(D,3,[],0.9);
% A = DMD.DMD.UX * DMD.DMD.Sig * DMD.DMD.
A = DMD.Xp * DMD.DMD.VX * pinv(DMD.DMD.Sig) * DMD.DMD.UX';
figDMDweigths(A);
% A = EDMD(D);
k = 2000; % k=round(0.0001*numel(A));
tic; [~, linIdx] = maxk(A(:), k); disp(toc);
E = zeros(k,2);
[E(:,1), E(:,2)] = ind2sub(size(A), linIdx);
GE = cell(k,2);
for i=1:k
    GE{i,1} = G(E(i,1));
    GE{i,2} = G(E(i,2));
end
figDMDgeneWeightHist(GE); % make a figure

% Load Hardwired Genome with full gene names
% HWG = load_HWG();                          % load HWG
load('C:\Users\picka\Documents\my_projects\DBTM\HardwiredGenome\Data\Processed\HWG\000\HWG.mat');
T = readtable('data/ensg2geneName.csv');   % load ensg to gene name map
ensg2name = containers.Map;                % build map so that HWG table may be populated faster
for i=1:height(T)
    e = string(T{i,1}{1});
    n = string(T{i,2}{1});
    ensg2name(e) = n;
end
HG = HWG.geneIndexTable;                   % populate HWG table with more gene names                    
for i=1:height(HG)
    if ensg2name.isKey(HG.("Stable ID")(i))
        n = ensg2name(HG.("Stable ID")(i));
        HG.("Gene Name")(i) = {n};
    end
end

% empty_count = sum(ismissing(HWG.geneIndexTable.("Gene Name")))
% empty_count = sum(ismissing(HG.("Gene Name")))

% sum(cellfun('isempty', HWG.geneIndexTable.("Gene Name")))
% sum(cellfun('isempty', HG.("Gene Name")))

% HWGUG = unique(HWG.geneIndexTable.("Gene Name"));
% HWGUG = unique(HG.("Gene Name"));

name2idx = containers.Map;                % build map so that HWG table may be populated faster
for i=1:height(HG)
    n = string(HG.("Gene Name"){i});
    if ~isempty(n)
        name2idx(n) = i;
    end
end

% Map MYOD data to the HWG
HWGIDXS = zeros(k,2);
m = 0;
for i=1:height(GE)
    n1 = string(GE{i,1}{1});
    n2 = string(GE{i,2}{1});
    if name2idx.isKey(n1) + name2idx.isKey(n2) == 2
        i1 = name2idx(n1);
        i2 = name2idx(n2);
        HWGIDXS(i,1) = i1;
        HWGIDXS(i,2) = i2;
    else
        HWGIDXS(i,1) = -1;
        HWGIDXS(i,2) = -1;
        m = m + 1;
    end
end

% Compate MYOD sparsity with HWG
S = 0;
H = 0;
M = 0;
for i=1:size(HWGIDXS,1)
    i1 = HWGIDXS(i,1);
    i2 = HWGIDXS(i,2);
    if i1 < 0
        S = S + 1;
        continue
    elseif (HWG.A(i1,i2) ~= 0) + (HWG.A(i2,i1) ~= 0) > 0
        H = H + 1;
    else
        M = M + 1;
    end
end

disp("H: " + string(H));
disp("M: " + string(M));
disp("S: " + string(S));
disp("DMD in HWG: " + string(H/M));

%% Sparse Promoting DMD

dmdsp()

%% Old
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
