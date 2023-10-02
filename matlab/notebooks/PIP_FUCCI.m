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

% clearvars -except GC DC replicates

%% Reduce the data to only include the transcription factors of the bunch
% Load Hardwired Genome with full gene names
load('C:\Users\picka\Documents\my_projects\DBTM\HardwiredGenome\Data\Processed\HWG\000\HWG.mat');
%{
T = readtable('data/ensg2geneName.csv');   % load ensg to gene name map
ensg2name = containers.Map;                % build map so that HWG table may be populated faster
for i=1:height(T)
    e = string(T{i,1}{1});
    n = string(T{i,2}{1});
    ensg2name(e) = n;
end
%}
HG = HWG.geneIndexTable;                   % populate HWG table with more gene names                    
hgname2idx = containers.Map;
for i=1:height(HG)
    if ~isempty(HG.("Gene Name"){i})
        hgname2idx(HG.("Gene Name"){i}) = i;
    end
end
TFS = zeros(numel(GC),1);
for i=1:size(GC)
    TFS(i) = HG.("Transcription Factor")(hgname2idx(GC{i}));
end
GT = GC(find(TFS ~= 0));
DT = DC(find(TFS ~= 0),:);

[s, GV] = hasnain2023(DT, replicates);
GC(s(1:5))

% sum(TFS)
% GC(find(TFS~=0))
% sum(HG.("Transcription Factor"))

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
[U,sigma,V] = svd(O');

figure; scatter(1:min(size(sigma)), diag(sigma),'.'); set(gca,'YScale','log');
title('Observability Matrix Scree Plot'); ylabel('Singular Value'); xlabel('Component');
figure; imagesc(A); title('DMD Linear Operator');
figure; plot(1:size(DC,2),DC(sensorIdxs,:));
title("PIP-FUCCI Sensors"); xlabel("Time"); ylabel("Expression"); legend(sensors);
figure; scatter(1:numel(A),sort(A(:),'descend'),'.'); title('DMD Weights');
xlabel('Component'); ylabel('Value'); % set(gca,'YScale','log');

%% Perform Sensor Selection by Hasnain and friends Nat. Comm. 2023
%
%   Note: GEM, PCNA, and CDT1 are never selected and the system does not
%   become observable according to this
[s, GV] = hasnain2023(DC, replicates);
dmd = shiftedDMD(DC,replicates,[],0.9);
A = dmd.Xp*dmd.DMD.VX*pinv(dmd.DMD.Sig)*dmd.DMD.UX';
for numSensors=1:numel(s)
    C = getC(s(1:numSensors),size(DC,1));
    O = obsv(A,C);
    r = rank(O);
    disp("# Sensors: " + string(length(unique(s(1:numSensors)))) + ", rank: " + string(r));
end
% GC(s(1:5))

%% What is the largest set of genes that can be observed from PIP-FUCCI
%
% Sunday Morning 10/1: 34 genes are fully observable: Can we do better/find
% a larger set? - JP

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

% perform DMD to construct A matrix
dmd = shiftedDMD(DC,3,[],0.9);
A = dmd.Xp*dmd.DMD.VX*pinv(dmd.DMD.Sig)*dmd.DMD.UX';

% visualize weights of A matrix
figure; scatter(1:numel(A),sort(A(:),'descend'),'.'); title('DMD Weights');
xlabel('Component'); ylabel('Value'); % set(gca,'YScale','log');

% Construct GRN from DMD A
thresh = 0.25;
Ag = (A > thresh) + (A < -1 * thresh); Ag = Ag + Ag';
g = graph(Ag);
figure; P = plot(g); P.NodeLabel = GC;

% Rank all genes according to their distances from CDT1, PCNA, and GEM
d = distances(g);
dists = zeros(size(DC,1), 3);
for s=1:3
    for i=1:size(DC,1)
        dists(i,s) = d(i,sensorIdxs(s));
    end
end
totalDist = sum(dists,2);
[~, idx] = sort(totalDist,'ascend');

% while the system is fully observable, add new genes until it is not fully
% observable
reducedSystem = sensorIdxs;
reducedSensors = [1 2 3];
r = 3;
i = 1;
while true
    DR = DC(reducedSystem,:);
    dmd = shiftedDMD(DR,3,[],0.9);
    Ar = dmd.Xp*dmd.DMD.VX*pinv(dmd.DMD.Sig)*dmd.DMD.UX';
    Cr = getC(reducedSensors, size(Ar,1));
    Or = obsv(Ar,Cr);
    r = rank(Or);
    if length(reducedSystem) == r
        % add a new vertex
        while find(reducedSystem == idx(i))
            i = i + 1;
        end
        reducedSystem = [reducedSystem idx(i)];
    else
        reducedSystem(end) = [];
        break
    end
end
disp(reducedSystem)
disp(GC(reducedSystem))

%% Which genes in GR are transcription factors
load('C:\Users\picka\Documents\my_projects\DBTM\HardwiredGenome\Data\Processed\HWG\000\HWG.mat');
HG = HWG.geneIndexTable;                   % populate HWG table with more gene names                    
name2idx = containers.Map;                % build map so that HWG table may be populated faster
for i=1:height(HG)
    n = string(HG.("Gene Name"){i});
    if ~isempty(n)
        name2idx(n) = i;
    end
end
TFs = HG.("Transcription Factor");
GRtf = zeros(length(GR),1);
for i=1:length(GR)
    idx = name2idx(GR{i});
    GRtf(i) = TFs(idx);
end

%% Make a graphic for the reduced system
DR = DC(reducedSystem,:);
dmd = shiftedDMD(DR,3,[],0.9);
A = dmd.Xp*dmd.DMD.VX*pinv(dmd.DMD.Sig)*dmd.DMD.UX';
GR = GC(reducedSystem);

% Set the figure properties
fig = figure;
set(fig, 'Units', 'inches', 'Position', [1, 1, 13, 6]); % Adjust the figure size

tiledlayout(1, 2);

% Create a heatmap
nexttile;
imagesc(A);
colormap('jet');
colorbar;
% set(gca, 'XTick', [], 'YTick', 1:34, 'YTickLabel', GR);
set(gca, 'XTick', 1:34, 'XTickLabel', GR, 'YTick', 1:34, 'YTickLabel', GR);
xtickangle(60); % Rotate X-axis labels for readability
title('Linearized Dynamics (DMD)', 'Interpreter', 'latex', 'FontSize', 14); % Use LaTeX for title

% Construct GRN from DMD A
nexttile;
% figure;
thresh = 2;
Ag = (A > thresh) + (A < -1 * thresh); Ag = Ag + Ag';
g = graph(Ag);

% Define which nodes to label and their categories
nodes_to_label = [1, 2, 3, find(GRtf == 1)'];
node_categories = {'Gene', 'Sensor', 'Transcription Factor'}; % Assign categories to nodes
node_colors = {'g', 'g', 'g', 'r', 'r'}; % Assign colors to categories

P = plot(g, 'NodeLabel', {}, 'MarkerSize', 10, 'LineWidth', 1.2); % Customize node appearance
% P.NodeLabel = {}; % Remove node labels for all nodes
P.NodeLabel = GR;

% Label nodes with categories and colors
for i=1:length(GR)
    if find(nodes_to_label == i)
        if i <= 3; color = 'g';
        else; color = 'r'; end
        highlight(P, i, 'NodeColor', color);
    else
        disp(GR(i))
        % label = "X";
        labelnode(P,[i],"")
        % P.NodeLabel{i} = label;
    end
end

title('Observable Cell Cycle Network', 'Interpreter', 'latex', 'FontSize', 14); % Use LaTeX for title

sgtitle('PIP-FUCCI Observability', 'Interpreter', 'latex', 'FontSize', 16); % Use LaTeX for super title
saveas(fig,'PIP-FUCCI.png')
%%
for i = 1:length(nodes_to_label)
    node = nodes_to_label(i);
    label = GR{node};
    % category = node_categories{i};
    color = node_colors{i};
    
    % P.NodeLabel(node) = GR(node);
    highlight(P, node, 'NodeColor', color);
end


title('Observable Cell Cycle Network', 'Interpreter', 'latex', 'FontSize', 14); % Use LaTeX for title

sgtitle('PIP-FUCCI Observability', 'Interpreter', 'latex', 'FontSize', 16); % Use LaTeX for super title

%%

% Create a legend
legend_categories = unique(node_categories);
legend_labels = cell(1, numel(legend_categories));
for i = 1:numel(legend_categories)
    legend_labels{i} = [legend_categories{i} ' (' node_colors{strcmp(node_categories, legend_categories{i})} ')'];
end
legend(legend_labels, 'Location', 'NorthEast');

% Create a legend manually
legend_categories = unique(node_categories);
legend_labels = cell(1, numel(legend_categories));

for i = 1:numel(legend_categories)
    category = legend_categories{i};
    color = node_colors{strcmp(node_categories, category)};
    legend_labels{i} = [category ' (' color ')'];
end

legend(legend_labels, 'Location', 'NorthEast');


title('Observable Cell Cycle Network', 'Interpreter', 'latex', 'FontSize', 14); % Use LaTeX for title

sgtitle('PIP-FUCCI Observability', 'Interpreter', 'latex', 'FontSize', 16); % Use LaTeX for super title


%%
% Visualize the weights
nexttile; scatter(1:numel(A),sort(A(:),'descend'),'.'); title('DMD Weights');
xlabel('Component'); ylabel('Value'); % set(gca,'YScale','log');

%% Correlations in the data
C = corr(DC');

figure; plot(sort(C(:),'descend'))
A = (C > 0.8); A = A - diag(diag(A));
figure; P = plot(graph(A)); P.NodeLabel = GC;

%% Perform DMD
dmd = shiftedDMD(DC,3,[],0.9);



