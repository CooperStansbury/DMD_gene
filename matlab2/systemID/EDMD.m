function [A] = EDMD(DataMatrix)
%EDMD Exact Dynamic Mode Decomposition
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 4, 2023

% Construct data matrices note, time snapshots are columns
    X  = DataMatrix(:,1:end-1);  
    Xp = DataMatrix(:,2:end);   

% Morse-Penrose Inverse
    Xi = pinv(X);

% Exact A
    A = Xp * Xi;


end

