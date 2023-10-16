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

filename = "hasnain2023_data_fc_gene_list.csv";
genes = readtable(filename,'Delimiter','\t');

filename = "hasnain2023_data_fc.csv";
T = readtable(filename);
D = T{2:end,2:end};

end

