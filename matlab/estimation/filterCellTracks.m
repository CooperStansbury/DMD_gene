function [cellTracks] = filterCellTracks(cellTracks, minSamples)
%FILTERCELLTRACKS Cleans the data to remove crazy/noisy/unwanted samples
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

%% Remove cells that have too few observations
ids = cellTracks.ID;
uidxs = unique(ids);
keep = [];
for i=1:length(uidxs)
    id = uidxs(i);
    idIdxs = find(ids==id);
    if numel(idIdxs) > minSamples
        keep = [keep; idIdxs];
    end
end

cellTracks = cellTracks(keep,:);

end

