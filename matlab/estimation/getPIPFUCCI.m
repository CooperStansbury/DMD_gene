function [outputs] = getPIPFUCCI(cellTracks, ids)
%GETPIPFUCCI Returns PIP-FUCCI signals from the cellt tracks table
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

if nargin == 1
    ids = unique(cellTracks.ID);
    % ids = mode(cellTracks.ID);
end

outputs = cell(length(ids),1);
systemOutput = [cellTracks.c0_intensity_mean cellTracks.c1_intensity_mean cellTracks.c2_intensity_mean];
for i=1:length(ids)
    id = ids(i);
    idIdxs = find(cellTracks.ID == id);
    Y = systemOutput(idIdxs,:);
    outputs{i} = Y;
end

end


