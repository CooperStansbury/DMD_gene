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

%% TO REMOVE
% Remove cells that ...
%    1. go backwards in cc
%    2. stall in cc
%    3. have c0 and c2 never turn on (the nucleus was not stained)
%    4. movement way to fast
%    5. based on distribution of nuclear area
%    6. position in the wound

%% Impute dummy variable
%   have a system to replace NAN values

end

