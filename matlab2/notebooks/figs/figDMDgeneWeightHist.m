function figDMDgeneWeightHist(GE)
%FIGDMDGENEWEIGHTHIST Summary of this function goes here
%   Detailed explanation goes here

%% Make plot
topGenes = 20;

GEE = GE;
for i=1:size(GEE,1)
    GEE{i,1} = (GEE{i,1}{1});
    GEE{i,2} = (GEE{i,2}{1});
end

GE1 = GEE(:,1);
GE2 = GEE(:,2);

% Count the frequency of each gene name
geneCounts = tabulate(GE1);

% Sort the geneCounts by frequency in descending order
geneCounts = sortrows(geneCounts, -2);

% Extract the top K gene names and their frequencies
topKGenes = geneCounts(1:topGenes, :);

% Extract gene names and their frequencies
geneNames1 = topKGenes(:, 1);
geneFrequencies1 = cell2mat(topKGenes(:, 2));

% Count the frequency of each gene name
geneCounts = tabulate(GE2);

% Sort the geneCounts by frequency in descending order
geneCounts = sortrows(geneCounts, -2);

% Extract the top K gene names and their frequencies
topKGenes = geneCounts(1:topGenes, :);

% Extract gene names and their frequencies
geneNames2 = topKGenes(:, 1);
geneFrequencies2 = cell2mat(topKGenes(:, 2));

figure;
subplot(1,2,1);
bar(geneFrequencies1)
xticks(1:topGenes)
xticklabels(geneNames1)
% xlabel('Gene Names')
ylabel('Frequency')
title('Most Susceptible')

subplot(1,2,2);
bar(geneFrequencies2)
xticks(1:topGenes)
xticklabels(geneNames2)
% xlabel('Gene Names')
ylabel('Frequency')
title('Most Influential')
sgtitle('DMD Pairwise Interaction Weights');

%% End section

% % Create a bar plot
% barh(geneFrequencies, 'b')
% yticklabels(geneNames)
% xlabel('Frequency')
% ylabel('Gene Names')
% title('Top 20 Gene Name Frequencies')
% 
% geneCounts1 = countcats(categorical(GEE(:,1)));
% geneCounts2 = countcats(categorical(GEE(:,2)));
% 
% % Create histograms for each column
% figure;
% 
% % Histogram for the first column (column1)
% subplot(2, 1, 1);
% bar(geneCounts1);
% title('Histogram of Gene Names (Column 1)');
% xlabel('Gene Names');
% ylabel('Frequency');
% 
% % Histogram for the second column (column2)
% subplot(2, 1, 2);
% bar(geneCounts2);
% title('Histogram of Gene Names (Column 2)');
% xlabel('Gene Names');
% ylabel('Frequency');
% 
% % Adjust subplot spacing
% sgtitle('Gene Name Histograms for Two Columns');
% 
% 
% 
% categorical(GEE(:,1))
% 
% % Extract the two columns of gene names
% column1 = GE(:, 1);
% column2 = GE(:, 2);
% 
% for i=1:size(GE,1)
%     column1{i} = string(column1{i});
%     column2{i} = string(column2{i});
% end
% 
% C1 = [];
% C2 = [];
% for i=1:size(GE)
%     C1 = [C1; column1{i}];
%     C2 = [C2; column1{i}];
% end
% 
% geneCounts1 = countcats(categorical(C1));
% geneCounts2 = countcats(categorical(C2));
% 
% % Count the occurrences of each gene name in each column
% geneCounts1 = countcats(categorical(column1));
% geneCounts2 = countcats(categorical(column2));

end

