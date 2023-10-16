function [velocities] = getVelocity(cellTracks, ids)
%GETVELOCITY Returns the velocity of individual cells
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

ts = 10; % 10 minutes between pip-fucci images

if nargin == 1
    ids = unique(cellTracks.ID);
    % ids = mode(cellTracks.ID);
end

velocities = cell(length(ids),1);
cellPositions = [cellTracks.x cellTracks.y];
for i=1:length(ids)
    id = ids(i);
    idIdxs = find(cellTracks.ID == id);
    cellMeasurements = cellPositions(idIdxs,:);
    velocity = zeros(height(cellMeasurements),1);
    for t=2:numel(velocity)
        % velocity = distance / time
        velocity(t) = sqrt((cellMeasurements(t,1) - cellMeasurements(t-1,1))^2 + ...
            (cellMeasurements(t,2) - cellMeasurements(t-1,2))^2) / ts;
    end 
    velocities{i} = velocity;
end

end

