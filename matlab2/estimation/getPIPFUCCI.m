function [outputs] = getPIPFUCCI(cellTracks, ids, type, gene)
%GETPIPFUCCI Returns PIP-FUCCI signals from the cellt tracks table
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

if nargin == 1
    ids = unique(cellTracks.ID);
    type = "mean";
    gene = "GEM";
    % ids = mode(cellTracks.ID);
elseif nargin == 2
    type = "log";
    gene = "GEM";
    % ids = mode(cellTracks.ID);
end

outputs = cell(length(ids),1);
if strcmp(type,"mean")
    systemOutput = [cellTracks.c0_intensity_mean cellTracks.c1_intensity_mean cellTracks.c2_intensity_mean];
else
    systemOutput = [cellTracks.c0_logFoldChange_processed cellTracks.c1_logFoldChange_processed cellTracks.c2_logFoldChange_processed];
end

for i=1:length(ids)
    id = ids(i);
    idIdxs = find(cellTracks.ID == id);
    Y = systemOutput(idIdxs,:);
    outputs{i} = Y; % normalizeSignals(Y, gene);
end

end


