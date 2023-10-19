%% Figures
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 19, 2023

%% PIP-FUCCI Signals
clear all; close all; clc;
minSamples = 10;
cellTracks = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
cellTracks = filterCellTracks(cellTracks,minSamples, true);
% id = mode(cellTracks.ID)
idIdxs = find(cellTracks.ID == mode(cellTracks.ID));

figure; tiledlayout(1,2); fs1 = 18; fs2 = 16;
tspan = (1/3)*(0:size(D,1)-1);
nexttile;
D = [cellTracks.c0_intensity_mean(idIdxs) cellTracks.c1_intensity_mean(idIdxs) cellTracks.c2_intensity_mean(idIdxs)];
plot(tspan,D(:,1),'r'); hold on;
plot(tspan,D(:,2),'k');
plot(tspan,D(:,3),'g');
title('PIP-FUCCI Signal','Interpreter', 'latex','FontSize',fs1);
xlabel("Time (hours)",'Interpreter', 'latex','FontSize',fs2)
ylabel("Signal Intensity",'Interpreter', 'latex','FontSize',fs2);
% legend(["CDT1","PCNA","GEM"]);
nexttile;
D = [cellTracks.c0_logFoldChange_processed(idIdxs) cellTracks.c1_logFoldChange_processed(idIdxs) cellTracks.c2_logFoldChange_processed(idIdxs)];
plot(tspan,D(:,1),'r'); hold on;
plot(tspan,D(:,2),'k');
plot(tspan,D(:,3),'g');
title('PIP-FUCCI Signal','Interpreter', 'latex','FontSize',fs1);
xlabel("Time (hours)",'Interpreter', 'latex','FontSize',fs2)
ylabel("Signal Intensity (log)",'Interpreter', 'latex','FontSize',fs2);
legend(["CDT1","PCNA","GEM"]);

%% PIP-FUCCI DMD
clear all; close all; clc;
minSamples = 10;
cellTracks = readtable('data/C1.tracks.full.csv'); % Read in cell tracks
cellTracks = filterCellTracks(cellTracks,minSamples, false);

figure; tiledlayout(1,4); fs1 = 18; fs2 = 16;
nexttile;
signals = getPIPFUCCI(cellTracks,unique(cellTracks.ID),"mean");
out = nonuniformShiftedDMD(signals,[],1);
A = out.DMD.A_bar;
imagesc(A);
title("Mean Signal",'Interpreter', 'latex','FontSize',fs1);
yticks([1 2 3]);
yticklabels({"CDT1","PCNA","GEM"});%,'Interpreter', 'latex','FontSize',fs2);
xticks([1 2 3]);
xticklabels({"CDT1","PCNA","GEM"});%,'Interpreter', 'latex','FontSize',fs2);

nexttile;
signals = getPIPFUCCI(cellTracks,unique(cellTracks.ID),"log");
out = nonuniformShiftedDMD(signals,[],1);
A = out.DMD.A_bar;
imagesc(A);
title("Log Signal",'Interpreter', 'latex','FontSize',fs1);
set(gca,'xtick',[]); set(gca,'ytick',[]);

sensorGenes = ["CDT1","PCNA","GEM"];
targetGenes = [];
for ds=1:2
    nexttile;
    if ds == 1
        [D,G,reps] = load2015(true); % Load data set
    elseif ds == 2
        [D,G,reps] = loadMYOD(); % Load data set
    end
    gene2idx = containers.Map;                              % Map gene names to indices
    for i=1:numel(G); gene2idx(string(G{i})) = i; end
    geneIdxs = [];                                     % Identify cell cycle indices from the data
    for i=1:length(sensorGenes)
        if gene2idx.isKey(string(sensorGenes(i)))
            geneIdxs(end+1) = gene2idx(string(sensorGenes(i)));
        end
    end
    D = D(geneIdxs,:);
    out = shiftedDMD(D,reps,[],1);
    A = out.DMD.A_bar;
    imagesc(A);
    if ds == 1
        title("2015",'Interpreter', 'latex','FontSize',fs1);
    else
        title("2018",'Interpreter', 'latex','FontSize',fs1);
    end
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
end

