function [s, GV] = hasnain2023(D, varargin)
%HASNAIN2023 Sensor selection from Learning perturbation-inducible cell states
% from observability analysis of transcriptome dynamics
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 13, 2023

% Set default values
gramT = size(D, 2);
dmdThresh = 0.9;
replicates = 1;

% Parse varargin
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'gramT')
        gramT = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'dmdThresh')
        dmdThresh = varargin{i + 1};
    elseif strcmpi(varargin{i}, 'reps')
        replicates = varargin{i + 1};
    else
        error('Unrecognized parameter name: %s', varargin{i});
    end
end

%% Dynamic Mode Decomposition
out = shiftedDMD(D,replicates,[],dmdThresh);
Atilda = out.DMD.UX' * out.Xp * out.DMD.VX * inv(out.DMD.Sig);

%% Reduced Grammarian
Ai = Atilda;
G = out.DMD.UX' * D(:,1) * D(:,1)' * out.DMD.UX;
for i=1:gramT
    G = G + Ai * out.DMD.UX' * D(:,1) * D(:,1)'*out.DMD.UX*Ai';
    Ai = Ai * Atilda;
end

[V,~] = eig(G);         % Eigen decomposition of reduced Grammarian
GV = out.DMD.UX * V;    % Approximate eigenvectors of full Grammarian

%% Sensor Rankings
[~, s] = max(GV);       % Take the largest index in each eigenvector

end

