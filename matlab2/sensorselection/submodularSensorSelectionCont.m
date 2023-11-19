function [S] = submodularSensorSelectionCont(A, varargin)
%SUBMODULARSENSORSELECTION 
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 15, 2023

% Set default values
gramT = 1;
maxSensors = 2;
n = size(A, 1);
outPath = "";
subCriteria = 1;

A = logm(A) / 8; % Continuous A from discrete A

% Parse varargin
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'maxSensors')
        maxSensors = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'subCriteria')
        subCriteria = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'outPath')
        outPath = varargin{i + 1};
    else
        error('Unrecognized parameter name: %s', varargin{i});
    end
end

%% Submodular optimization
S = [];                     % selected sensors
R = 1:n;                    % remaining sensors

% while selecting more sensors
while numel(S) < maxSensors
    M = zeros(numel(R), 1);     % save scores for each sensor
    % try each of the remaining sensors
    parfor i=1:numel(R)
        disp(string(i) + "/" + string(numel(R)));
        vx = R(i);                  % pick an unused sensor
        C = getC([S vx], n);        % create C matrix

        G = lyap(A, C'*C);
        
        if subCriteria == 1
            M(i) = trace(G);      % Four measures of submodularity
        elseif subCriteria == 2
            M(i) = trace(inv(G));
        elseif subCriteria == 3
            M(i) = log(det(G));
        elseif subCriteria == 4
            M(i) = rank(G);
        end
    end
    [~, vx] = max(M); % find highest weighted next sensor
    S = [S R(vx)];    % select the next sensor
    R(vx) = [];       % remove sensor from remaining vertices

    % save results only if a outPath has been specified
    if ~strcmp(outPath, "")
        out = struct;
        out.subCriteria = subCriteria;
        out.M = M;
        out.S = S;
        out.R = R;
        out.n = n;
        filePath = outPath + "soss_itr" + string(length(S)) + "_criteria" + string(subCriteria) + "_gramT" + string(gramT) + ".mat";
        save(filePath, 'out', '-v7.3');
    end
end


end

