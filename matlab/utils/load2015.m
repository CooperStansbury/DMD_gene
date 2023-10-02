function [D,genes,reps] = load2015()
%LOADMYOD Loads 2015 data
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 2, 2023
T = readtable('chen2015normalized.txt');
genes = T.geneName;
D = T(:,2:end);
D = D{:,:};
reps = 2;

D1 = zeros(size(D,1),size(D,2)/2);                         % Create matrices to store data
D2 = zeros(size(D,1),size(D,2)/2);
for t=1:9                                          % Read each timepoint
    D1(:,t) = D(:,2*(t-1) + 1);
    D2(:,t) = D(:,2*(t-1) + 1);
end

D = [D1 D2];                                     % Organize the data

end

