%% Test Hasnain 2023
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 15, 2023

clear; close all; clc;

T = readtable('data/tpm_0.5_gene_list.csv');
[D, G] = loadRNAseq(1, T.gene_name);
S = hasnain2023(D, 'reps', 2, 'dmdThresh', 0.9, 'gramT', 20);
G(S)


%% Scratch
D
targetGenes = cellstr(T.gene_name)
geneNamesArray = table2cell(T.gene_name);
geneNamesArray = T.gene_name;
geneNamesArray = geneNamesArray(:);