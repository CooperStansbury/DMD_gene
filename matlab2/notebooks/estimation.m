%% Estimation
%
%   Problems:
%       - interpolate time series RNAseq data
%       - KF time series RNAseq data and show convergence
%       - interpolate PIP-FUCCI data (signal and position)
%       - KF PIP-FUCCI data and show convergence
%       - KF RNAseq data from PIP-FUCCI signals and show convergence
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 13, 2023

close all; clear; clc;

%% Interpolate PIP-FUCCI data from RNAseq
ds = 1;
sensorGenes = ["PCNA","CDT1","GEM"];
geneColor = ['k','r','g'];

[D, G] = loadRNAseq(ds, sensorGenes);
reps = 1; % size(D,3);
T = size(D,2);
n = size(D,1);

sfreq = 8;
delta = 1/3;
tts = sfreq * [0:T-1];
ttd = 0:delta:max(tts);


r=1;
    % Interpolate data forward
    DD = D(:,:,r);
    Ahat = EDMD(DD);      % Ahat at sfreq    
    AhatDelta = expm(logm(Ahat) * delta / sfreq);  % AhatDelta interpolates data at a higher frequency
    Di = [];
    for t=1:T-1
        Dt = DD(:, t);
        for i = 1:sfreq/delta-1
            Dt = [Dt AhatDelta * Dt(:,end)];
        end
        Di = [Di Dt];
    end
    Df = [Di DD(:,end)];

    % Interpolate data backward
    DD = flip(D(:,:,r), 2);
    Ahat = EDMD(DD);      % Ahat at sfreq    
    AhatDelta = expm(logm(Ahat) * delta / sfreq);  % AhatDelta interpolates data at a higher frequency
    Di = [];
    for t=1:T-1
        Dt = DD(:, t);
        for i = 1:sfreq/delta-1
            Dt = [Dt AhatDelta * Dt(:,end)];
        end
        Di = [Di Dt];
    end
    Db = [Di DD(:,end)];
    Db = flip(Db, 2);

figure;
tiledlayout(6,1);

nexttile;
for g=1:3
    scatter(tts, D(g,:,r), geneColor(g), 'o', 'filled'); hold on;
    scatter(ttd, Df(g,:),  geneColor(g), '_'); hold on;
end
title('Forward Interpolation'); ylabel('Expression');

nexttile;
for g=1:3
    scatter(tts, D(g,:,r), geneColor(g), 'o', 'filled'); hold on;
    scatter(ttd, Db(g,:),  geneColor(g), '|'); hold on;
end
title('Backward Interpolation'); ylabel('Expression');

nexttile;
for g=1:3
    scatter(tts, D(g,:,r), geneColor(g), 'o', 'filled'); hold on;
    scatter(ttd, Df(g,:),  geneColor(g), '_'); hold on;
    scatter(ttd, Db(g,:),  geneColor(g), '|'); hold on;
end
title('Forward-Backward Interpolation'); ylabel('Expression');

nexttile;
for g=1:3
    scatter(tts, D(g,:,r), geneColor(g), 'o', 'filled'); hold on;
    avgInt = (Df(g,:) + Db(g,:)) / 2;
    scatter(ttd, avgInt,  geneColor(g), '.'); hold on;
end
title('Unweighted Avg. Interpolation'); ylabel('Expression');

nexttile;
fvals = zeros(length(ttd), 1);
for g=1:3
    scatter(tts, D(g,:,r), geneColor(g), 'o', 'filled'); hold on;
    interp = zeros(numel(Df(g,:)), 1);
    for i=1:numel(interp)
        tival = mod(i, sfreq / delta); disp(tival);
        % If the data point is observed
        if tival == 1
            interp(i) = D(g, i / (sfreq / delta),r);
            fvals(i) = 0;
        else
            f = 1 - (tival / (sfreq/delta));
            interp(i) = ((f * Df(g,i)) + ((1-f) * Db(g,i)));
            if g==1; fvals(i) = f; end
        end
    end
    scatter(ttd, interp,  geneColor(g), '.'); hold on;
end
title('Time Weighted Interpolation'); ylabel('Expression');

nexttile;
scatter(ttd, fvals');

L1 = plot(nan, nan, 'color', 'k');
L2 = plot(nan, nan, 'color', 'r');
L3 = plot(nan, nan, 'color', 'g');
L4 = plot(nan, nan, '_', 'color', 'k');
L5 = plot(nan, nan, '|', 'color', 'k');
legend([L1, L2, L3, L4, L5], {'PCNA', 'CDT1', 'GEM', 'Forward Interpolation', 'Backward Interpolation'})
% title('Interpolating PIP-FUCCI Signals from RNAseq with FB DMD Models');
xlabel('Hours');


%% Interpolate PIP-FUCCI data from microscopy

sensorGenes = ["PCNA","CDT1","GEM"];
geneColor = ['k','r','g'];

cellTracks = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
cellTracks = filterCellTracks(cellTracks, 10);
signals = getPIPFUCCI(cellTracks,unique(cellTracks.ID),'log');

instance = 2;
sfreq = 1/3;
delta = 1/6;
figure; tiledlayout(2,1);

nexttile;
Y = signals{instance};
tt = (1/3) * (0:numel(Y(:,1))-1);
for g=1:3
    scatter(tt, Y(:,g), geneColor(g)); hold on;
end
xlabel('Hours');
ylabel('Expression');
title('PIP-FUCCI Microscopy Signals');

nexttile;
Y = Y';
T = size(Y,2);

tts = sfreq * [0:T-1];
ttd = 0:max(tts);

DD = Y;
Ahat = EDMD(DD);      % Ahat at sfreq
AhatDelta = expm(logm(Ahat) * delta / sfreq);  % AhatDelta interpolates data at a higher frequency
Di = [];
for t=1:T-1
    Dt = DD(:, t);
    for i = 1:sfreq/delta-1
        Dt = [Dt AhatDelta * Dt(:,end)];
    end
    Di = [Di Dt];
end
Df = [Di DD(:,end)];

% Interpolate data backward
DD = flip(Y, 2);
Ahat = EDMD(DD);      % Ahat at sfreq    
AhatDelta = expm(logm(Ahat) * delta / sfreq);  % AhatDelta interpolates data at a higher frequency
Di = [];
for t=1:T-1
    Dt = DD(:, t);
    for i = 1:sfreq/delta-1
        Dt = [Dt AhatDelta * Dt(:,end)];
    end
    Di = [Di Dt];
end
Db = [Di DD(:,end)];
Db = flip(Db, 2);

for g=1:3
    scatter(tts, Y(g,:), geneColor(g), 'o', 'filled'); hold on;
    scatter(ttd, Df(g,:),  geneColor(g), '_'); hold on;
    scatter(ttd, Db(g,:),  geneColor(g), '|'); hold on;
end

%% scratch
T = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
TT = [T.ID T.c0_intensity_mean T.c1_intensity_mean T.c2_intensity_mean];
% TT(TT(:,1) == 11,:) = [];
id = mode(TT(:,1));
TT(TT(:,1) ~= id,:) = [];

tt = (1/3) * (0:height(TT)-1);
figure;
for g=1:3
    scatter(tt, TT(:,g+1)); hold on;
end



%% scratch
if ds == 1
    [Dm,G,reps] = load2015(false); % Load data set
elseif ds == 2
    [Dm,G,reps] = loadMYOD(); % Load data set
end
T = size(Dm,2)/reps;
D = zeros(size(Dm,1), size(Dm,2)/reps, reps);
for i=1:reps
    D(:,:,i) = Dm(:,(i-1)*T+1:i*T);
end
% subset the data to only include sensor and target genes
gene2idx = containers.Map;                              % Map gene names to indices
for i=1:numel(G); gene2idx(string(G{i})) = i; end
geneIdxs = [];                                     % Identify cell cycle indices from the data
for i=1:length(sensorGenes)
    if gene2idx.isKey(string(sensorGenes(i)))
        geneIdxs(end+1) = gene2idx(string(sensorGenes(i)));
    end
end
Dwhole = D;
D = D(geneIdxs,:,:);
G = G(geneIdxs);

scatter3(signals{13}(:,1),signals{13}(:,2),signals{13}(:,3),'.');

