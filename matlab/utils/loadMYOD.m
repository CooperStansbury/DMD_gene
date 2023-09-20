function [D,genes] = loadMYOD()
%LOADMYOD Loads 2018 MYOD data
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
for t=2:16                                          % Read each timepoint
    D1(t,:) = T{3*(t-1) + 1, 2:size(T,2)};
    D2(t,:) = T{3*(t-1) + 2, 2:size(T,2)};
    D3(t,:) = T{3*(t-1) + 3, 2:size(T,2)};
end

D1 = D1'; D2 = D2'; D3 = D3';
D = [D1 D2 D3];                                     % Organize the data

end

