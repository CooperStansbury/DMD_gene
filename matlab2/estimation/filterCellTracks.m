function [cellTracks] = filterCellTracks(cellTracks, minSamples, isPlot)
%FILTERCELLTRACKS Cleans the data to remove crazy/noisy/unwanted samples
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

if nargin == 2
    isPlot = false;
end

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

%% Remove cells where one channel has near zero varriance
ids = cellTracks.ID;
uidxs = unique(ids);
SC0 = []; SC1 = []; SC2 = [];
keep = [];
for i=1:length(uidxs)
    id = uidxs(i);
    idIdxs = find(ids == id);
    c0 = cellTracks.c0_logFoldChange_processed(idIdxs);
    c1 = cellTracks.c1_logFoldChange_processed(idIdxs);
    c2 = cellTracks.c2_logFoldChange_processed(idIdxs);
    % Treshhold is based on a figure JP made. It can/should be modified
    if std(c0) > 0.1 && std(c1) > 0.1 & std(c2) > 0.1
        keep = [keep; idIdxs];
    end
    if isPlot; SC0 = [SC0; std(c0)]; SC1 = [SC1; std(c1)]; SC2 = [SC2; std(c2)]; end
end
cellTracks = cellTracks(keep,:);
if isPlot
    figure; histogram(SC0); hold on; histogram(SC1); histogram(SC2);
    legend(["PCNA","CDT1","GEM"]); title("Standard Deviation in PIP-FUCCI Log Fold Change Signals");
    xlabel("Standard Deviation"); ylabel("Unique Cell");
end

%% Remove cells with dummy variables (aan artifact of the image processing)
% dummy = cellTracks.dummy;

%% Impute dummy variable
%   have a system to replace NAN values

end

