%% Sensor Selection
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 17, 2023

clear; close all; clc

%% 2018 MYOD Data Set
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
clearvars -except D genes                           % Remove all other variables

% Perform sensor selection
[S, GV] = hasnain2023(D);
genes(S)

%% 2023 Hasnain Data set

filename = "hasnain2023_data_fc_gene_list.csv";
genes = readtable(filename,'Delimiter','\t');

filename = "hasnain2023_data_fc.csv";
T = readtable(filename);
D = T{2:end,2:end};
clearvars -except D genes                           % Remove all other variables

% Perform sensor selection
[S, GV] = hasnain2023(D);
genes.GeneName(S)
