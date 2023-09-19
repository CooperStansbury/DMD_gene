function [D,genes] = loadHasnain2023()
%LOADHASNAIN2023 Loads data from Hasnain 2023
%
% Reference: Learning Perturbation-Induced Cell States from Observability 
%            Analysis of Transcriptome Dynamics
%
% Outputs:
%     - D: Data matrix of gene expression x time. There are 2 replicates 
%          concatenated horizonatally where time point 1 for rep. 1 is column
%          1 and time point 1 for rep. 2 is column 10
%     - genes: an array containing gene names corresponding to each row in D
%
% Data Processing: These data are obtained directly from the code published
%                  with the paper above using the same preprocessing steps. In
%                  particular, the first few cells in main.ipynb were run
%                  until data_fc and genes_keep were obtained. These objects
%                  were written to the corresponding .csv files loaded here.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 18, 2023

T = readtable("time_series_rna_fold_changes.csv");  % Read .csv file
fid = fopen("time_series_rna_fold_changes.csv",'r');
genes = fgetl(fid); genes = split(genes, ',');
genes = genes(2:end);
D1 = zeros(15,size(T,2)-1);                         % Create matrices to store data
D2 = zeros(15,size(T,2)-1);
D3 = zeros(15,size(T,2)-1);
for t=1:15                                          % Read each timepoint
    D1(t,:) = T{3*(t-1) + 1, 2:size(T,2)};
    D2(t,:) = T{3*(t-1) + 2, 2:size(T,2)};
    D3(t,:) = T{3*(t-1) + 3, 2:size(T,2)};
end

D1 = D1'; D2 = D2'; D3 = D3';
D = [D1 D2 D3];                                     % Organize the data

end

