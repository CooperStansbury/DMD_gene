function [Y] = normalizeSignals(Y,gene)
%NORMALIZESIGNALS 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 9, 2023

% gene = "GEM";
switch gene
    case "CDT1"
        idx = 1;
    case "PCNA"
        idx = 2;
    case "GEM"
        idx = 3;
end

[~, t] = max(Y(:,idx));
timeMax = Y(t,:);
Y = Y ./ timeMax;

end

