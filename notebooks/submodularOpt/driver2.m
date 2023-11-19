function driver2(Apath, ms, op, sc)
%DRIVER FOR SUBMODULAR OPTIMIZATION WITH CONTINUOUS RELAXATION
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 16, 2023

disp('Job Starting');
disp(ms);
disp(op);
disp(sc);
disp(Apath);

% Add paths
addpath(genpath('/home/jpic/DMD_gene/matlab2/'))

% Load data
load(Apath)

% Call optimizer
submodularSensorSelectionCont(A, 'maxSensors', ms, 'outPath', op,'subCriteria', sc)

end