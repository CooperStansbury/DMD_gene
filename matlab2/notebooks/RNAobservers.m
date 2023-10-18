%% RNA Observers
%
%   This file uses DMD + LQE/Kalman estimation to estimate gene expression
%   from a subset of genes
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 17, 2023

clear; close all; clc;
isPlot = 1;

ds = 2;
sensorGenes = ["PCNA","CDT1","GEM"];
% nSensors = 10;
% targetGenes = [];
% Ovservable genes from CDT1, PCNA, and GEM according to PIP-FUCCI.m
targetGenes = ["PPP2R5B","CCND1","CCND2","CCND3","CDK4","CDK6","RB1","RBL1","RBL2","ABL1","HDAC1","HDAC2","E2F1","E2F2","E2F3","E2F4","E2F5"];%,"TFDP1","TFDP2","GSK3B","TGFB1","TGFB3","SMAD2","SMAD3","SMAD4"];

if ds == 1
    [D,G,reps] = load2015(false); % Load data set
    % sensorGenes = ["PSRC1", "UBE2C", "FAM72A", "KIF20A", "FAM72C", "GAS2L3", "CENPA", "PLK1", "TOP2A", "CDCA3", "ASPM", "DEPDC1", "FAM72B", "AURKA", "KIF23", "KIF14", "FAM72D", "ERN2", "CDC25C", "ANLN", "NDC80", "GTSE1", "GNB3", "HJURP", "IQGAP3", "CDCA2", "CDK1", "CKAP2L", "NEK2", "DLGAP5", "CDCA8", "FBXO43", "KIF2C", "TROAP", "UNC13A", "KPNA2", "ARHGAP11A", "NUF2", "TTK", "SGO2", "RACGAP1", "HMMR", "AURKB", "TPX2", "INCENP", "CCNB1", "DEPDC1B", "DHRS2", "BUB1", "SPC25", "ECT2", "CCNF", "SKA1", "BIRC5", "KIF11", "CDC20", "ARHGAP11B", "KIF18B", "KNSTRN", "LCE1E", "CCNA2", "SPAG5", "KIFC1", "CCNB2", "MYLK2", "KIF4A", "CEP70", "KIF4B", "ARHGAP19", "PIMREG", "CENPI", "KNL1", "PRR11", "PBK", "CENPE", "NCAPG", "BUB1B", "ESPL1", "CEP55", "KIF18A", "DIAPH3", "TACC3", "DRP2", "NT5DC4", "STXBP5L", "CCDC150", "ARL6IP1", "SORBS1", "APOBEC3B", "KIF15", "H2BC18", "CKS2", "CIT", "PRC1", "HMGB2", "SGO1", "CDKN3", "CEP128", "DNAJC22", "SPC24", "NOSTRIN", "NEIL3", "LIPI", "PIF1", "NUSAP1", "CA9", "CKAP2", "CKS1B", "H3C8", "PTTG2", "TUBA1C", "FOXM1", "NCAPH", "SFN", "PRR5L", "MYORG", "MELK", "UBASH3A", "H3C3", "MKI67", "KIF22", "GPR19", "PADI2", "CATSPERE", "THSD4", "ARHGAP19-SLIT1", "KIF20B", "SDSL", "H3C2", "SRGAP2", "CYP4F3", "NET1", "NOTCH2NLC", "CIP2A", "BORA", "ATP8B3", "SPDL1", "STIL", "SKA3", "KCNJ13", "C12orf71", "KIAA1671", "PLK4", "SAPCD2", "H2AC20", "SHCBP1", "MTFR2", "TLCD3B", "MID1", "CKAP5", "CENPF", "FAM83D", "KCTD14", "ALAS2", "ACKR2", "MAMDC2", "RASGRP4", "EPHA3", "BRD8", "SHB", "PARPBP", "NCAPD2", "TRAIP", "ERCC6L", "PMCH", "RNF26", "H1-5", "MYADML2", "OIP5", "ARHGEF39", "MYO1H", "USP2", "ALOX12B", "C7orf25", "KIF7", "C18orf54", "IRF6", "C2orf78", "UBE2S", "GDA", "NECTIN4", "POLQ", "PADI4", "H2BC12", "TMPRSS5", "G2E3", "INSL6", "SMC4", "TICRR", "RDM1", "C21orf58", "H2BC7", "REXO5", "ZNF850", "NCAPG2", "POC1A", "CATSPER3", "DSG1", "SDC4", "TRIM59", "GNAT3", "CENPW", "CLGN", "SRGAP2C", "LRRC49", "OR56B1", "HYLS1", "TAS2R50", "TUBB", "PTPN14", "ANG", "DNMT3B", "RERGL", "H2AX", "TEX49", "MAD2L1", "DBF4B", "HASPIN", "FRMPD2", "GGH", "CCDC163", "PTTG1", "ALOX15", "TRIP13", "HSPA1L", "HROB", "MPP4", "ARVCF", "FANCD2", "HTRA4", "PNCK", "SH3PXD2B", "FXYD6", "CDC25B", "NLRP9", "DSCAML1", "H3C7", "OR2B6", "PHYHIP", "HMGB3", "STEAP2", "CDKN2C", "H2BC11", "TMEM45A", "CYP4F2", "TUBB4B", "LCE1F", "ARHGAP29", "COL4A2", "CDCA5", "SLC52A1", "DIXDC1", "DLG3", "POTEI", "PAPPA", "RRM2", "POLH", "KCNV2", "CDK15", "LRRC3", "DNAH6", "PRELID2", "FRY", "FANCD2OS", "CDKN2D", "C5orf34", "KCNA6", "MZT1", "SHTN1", "SRGAP2B", "ZWINT", "PDE8B", "PRCD", "GXYLT2", "SMTN", "CCDC15", "TMEM88B", "TNS1", "SERPINE1", "ZNF165", "ODAD3", "METTL4", "KCNJ15", "OGN", "CSPG4", "CASP14", "REEP4", "AXDND1", "H4C6", "DBF4", "FAM166A", "NKAPL", "ADAMTSL1", "SUV39H1", "GPSM2", "SCG5", "MORN3", "STMN1", "TUBA1B", "ESCO2"];
    % sensorGenes = sensorGenes(1:nSensors);
elseif ds == 2
    [D,G,reps] = loadMYOD(); % Load data set
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
for i=1:length(targetGenes)
    if gene2idx.isKey(string(targetGenes(i)))
        geneIdxs(end+1) = gene2idx(string(targetGenes(i)));
    end
end
Dwhole = D;
D = D(geneIdxs,:);
G = G(geneIdxs);

% Construct system matrices (A,C)
out = shiftedDMD(D, reps, [], 1);
A = out.DMD.A_bar;
C = getC(1:3, length(geneIdxs)); C = full(C);

if isPlot
    E = eig(A); r = real(E); i = imag(E); theta = linspace(0, 2*pi, 1000); x = cos(theta); y = sin(theta);
    figure; plot(x, y); hold on; scatter(r,i); xlabel('Real ($\lambda$)','Interpreter', 'latex');
    ylabel('Imaginary ($\lambda$)','Interpreter', 'latex'); title('DMD Eigenvalues','Interpreter', 'latex');
end

% Observability over the 3 replicates from 2018 MYOD
itrs = 1000;

B0 = zeros(size(A,1),1);    % Construct kalman filter
D0 = zeros(size(C,1),1);
sym = ss(A,B0,C,D0);
[kalmf,~,~] = kalman(sym,1e-4,1e-4);

nSensors = numel(sensorGenes);
T = size(D,2) / reps;
tspan = 1:T-1;
E = cell(3,1);
for r=1:reps
    xtrue = D(:,T*(r-1)+2:T*r);
    ytrue = xtrue(1:nSensors,:);
    ER = [];
    for i=1:itrs
        x0 = xtrue(:,1); % Set initial entry
        ind_1 = randi([1 size(Dwhole,1)],[numel(x0) - (nSensors), 1]);
        x0(nSensors + 1:end) = Dwhole(ind_1,2);
        % x0(nSensors + 1:end,:) = rand(numel(x0) - nSensors,1);
        xhat = lsim(kalmf,ytrue,tspan,x0); xhat = xhat';
        xhat = xhat(length(sensorGenes) + 1:end,:);
        ei = [];
        for t=1:size(xtrue,2)
            xh = xhat(:,t) / sum(xhat(:,t));
            xt = xtrue(:,t) / sum(xtrue(:,t));
            e = norm(xh - xt);   % we could use KL
            ei = [ei e];
        end
        ER = [ER; ei];
    end
    E{r} = ER;
end

% Random guessing
ER = [];
for i=1:itrs
    ei = [];
    for t=1:size(xtrue,2)
        % xh = xhat(:,t) / sum(xhat(:,t));
        xh = rand(); xh = xh / sum(xh);
        xt = xtrue(:,t) / sum(xtrue(:,t));
        e = norm(xh - xt);   % we could use KL
        ei = [ei e];
    end
    ER = [ER; ei];
end
E{r+1} = ER;

% Visualize the errors
figure;
tspan = 8*(1:size(ER,2));
for r=1:reps+1
    ER = E{r};
    m = mean(ER);
    err = std(ER) / sqrt(itrs);
    errorbar(tspan,m,err);
    if r == 1; hold on; end
end
title('LQE/Kalman Filter Estimator Error','Interpreter', 'latex');
legend(["Replicate 1", "Replicate 2", "Replicate 3","Random"]);
xlabel("Time (hours)",'Interpreter', 'latex')
ylabel("Estimate Error",'Interpreter', 'latex');

%% X0 estimation from observability matrix

clear; close all; clc;
isPlot = 1;

ds = 2;
% Ovservable genes from CDT1, PCNA, and GEM according to PIP-FUCCI.m
sensorGenes = ["PCNA","CDT1","GEM"]; nSensors = numel(sensorGenes);
% nSensors = 10;
% targetGenes = [];
targetGenes = ["PPP2R5B","CCND1","CCND2","CCND3","CDK4","CDK6","RB1","RBL1","RBL2","ABL1","HDAC1","HDAC2","E2F1","E2F2","E2F3","E2F4","E2F5"];%,"TFDP1","TFDP2","GSK3B","TGFB1","TGFB3","SMAD2","SMAD3","SMAD4"];

if ds == 1
    [D,G,reps] = load2015(false); % Load data set
    % sensorGenes = ["PSRC1", "UBE2C", "FAM72A", "KIF20A", "FAM72C", "GAS2L3", "CENPA", "PLK1", "TOP2A", "CDCA3", "ASPM", "DEPDC1", "FAM72B", "AURKA", "KIF23", "KIF14", "FAM72D", "ERN2", "CDC25C", "ANLN", "NDC80", "GTSE1", "GNB3", "HJURP", "IQGAP3", "CDCA2", "CDK1", "CKAP2L", "NEK2", "DLGAP5", "CDCA8", "FBXO43", "KIF2C", "TROAP", "UNC13A", "KPNA2", "ARHGAP11A", "NUF2", "TTK", "SGO2", "RACGAP1", "HMMR", "AURKB", "TPX2", "INCENP", "CCNB1", "DEPDC1B", "DHRS2", "BUB1", "SPC25", "ECT2", "CCNF", "SKA1", "BIRC5", "KIF11", "CDC20", "ARHGAP11B", "KIF18B", "KNSTRN", "LCE1E", "CCNA2", "SPAG5", "KIFC1", "CCNB2", "MYLK2", "KIF4A", "CEP70", "KIF4B", "ARHGAP19", "PIMREG", "CENPI", "KNL1", "PRR11", "PBK", "CENPE", "NCAPG", "BUB1B", "ESPL1", "CEP55", "KIF18A", "DIAPH3", "TACC3", "DRP2", "NT5DC4", "STXBP5L", "CCDC150", "ARL6IP1", "SORBS1", "APOBEC3B", "KIF15", "H2BC18", "CKS2", "CIT", "PRC1", "HMGB2", "SGO1", "CDKN3", "CEP128", "DNAJC22", "SPC24", "NOSTRIN", "NEIL3", "LIPI", "PIF1", "NUSAP1", "CA9", "CKAP2", "CKS1B", "H3C8", "PTTG2", "TUBA1C", "FOXM1", "NCAPH", "SFN", "PRR5L", "MYORG", "MELK", "UBASH3A", "H3C3", "MKI67", "KIF22", "GPR19", "PADI2", "CATSPERE", "THSD4", "ARHGAP19-SLIT1", "KIF20B", "SDSL", "H3C2", "SRGAP2", "CYP4F3", "NET1", "NOTCH2NLC", "CIP2A", "BORA", "ATP8B3", "SPDL1", "STIL", "SKA3", "KCNJ13", "C12orf71", "KIAA1671", "PLK4", "SAPCD2", "H2AC20", "SHCBP1", "MTFR2", "TLCD3B", "MID1", "CKAP5", "CENPF", "FAM83D", "KCTD14", "ALAS2", "ACKR2", "MAMDC2", "RASGRP4", "EPHA3", "BRD8", "SHB", "PARPBP", "NCAPD2", "TRAIP", "ERCC6L", "PMCH", "RNF26", "H1-5", "MYADML2", "OIP5", "ARHGEF39", "MYO1H", "USP2", "ALOX12B", "C7orf25", "KIF7", "C18orf54", "IRF6", "C2orf78", "UBE2S", "GDA", "NECTIN4", "POLQ", "PADI4", "H2BC12", "TMPRSS5", "G2E3", "INSL6", "SMC4", "TICRR", "RDM1", "C21orf58", "H2BC7", "REXO5", "ZNF850", "NCAPG2", "POC1A", "CATSPER3", "DSG1", "SDC4", "TRIM59", "GNAT3", "CENPW", "CLGN", "SRGAP2C", "LRRC49", "OR56B1", "HYLS1", "TAS2R50", "TUBB", "PTPN14", "ANG", "DNMT3B", "RERGL", "H2AX", "TEX49", "MAD2L1", "DBF4B", "HASPIN", "FRMPD2", "GGH", "CCDC163", "PTTG1", "ALOX15", "TRIP13", "HSPA1L", "HROB", "MPP4", "ARVCF", "FANCD2", "HTRA4", "PNCK", "SH3PXD2B", "FXYD6", "CDC25B", "NLRP9", "DSCAML1", "H3C7", "OR2B6", "PHYHIP", "HMGB3", "STEAP2", "CDKN2C", "H2BC11", "TMEM45A", "CYP4F2", "TUBB4B", "LCE1F", "ARHGAP29", "COL4A2", "CDCA5", "SLC52A1", "DIXDC1", "DLG3", "POTEI", "PAPPA", "RRM2", "POLH", "KCNV2", "CDK15", "LRRC3", "DNAH6", "PRELID2", "FRY", "FANCD2OS", "CDKN2D", "C5orf34", "KCNA6", "MZT1", "SHTN1", "SRGAP2B", "ZWINT", "PDE8B", "PRCD", "GXYLT2", "SMTN", "CCDC15", "TMEM88B", "TNS1", "SERPINE1", "ZNF165", "ODAD3", "METTL4", "KCNJ15", "OGN", "CSPG4", "CASP14", "REEP4", "AXDND1", "H4C6", "DBF4", "FAM166A", "NKAPL", "ADAMTSL1", "SUV39H1", "GPSM2", "SCG5", "MORN3", "STMN1", "TUBA1B", "ESCO2"];
    % sensorGenes = sensorGenes(1:nSensors);
elseif ds == 2
    [D,G,reps] = loadMYOD(); % Load data set
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
for i=1:length(targetGenes)
    if gene2idx.isKey(string(targetGenes(i)))
        geneIdxs(end+1) = gene2idx(string(targetGenes(i)));
    end
end
Dwhole = D;
D = D(geneIdxs,:);
G = G(geneIdxs);

% Construct system matrices (A,C)
out = shiftedDMD(D, reps, [], 1);
A = out.DMD.A_bar;
C = getC(1:3, length(geneIdxs)); C = full(C);

% Set up least squares problem
O = obsvt(A,C,)
O = obsv(A,C);
Y = D(1:nSensors,:);
y = reshape(Y,[numel(Y) 1]);

%% state space model

B0 = zeros(size(A,1),1);    % Construct kalman filter
D0 = zeros(size(C,1),1);
sym = ss(A,B0,C,D0);
[kalmf,~,~] = kalman(sym,1e-4,1e-4);

% kalmf.OutputGroup
% kalmf.A

% Apply the kalman filter
tspan = 1:1:15;
x0 = D(:,2); x0(length(sensorGenes) + 1:end,:) = rand(numel(x0) - numel(sensorGenes),1);
y = lsim(kalmf,zeros(numel(tspan),3),tspan,x0); y = y';
xhat = y(length(sensorGenes) + 1:end,:);
xtrue = D(:,2:16);

% Error analysis
E = [];
for t=1:size(xtrue)
    e = norm(xhat(:,t) - xtrue(:,t));
    E = [E e];
end
figure; plot(E);


%%
[kalmf,L,P] = kalman(sym,0,0)
kalman(sym)
% Q = cov(D'); R = cov(D(1:length(sensorGenes),:)');