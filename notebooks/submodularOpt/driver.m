function driver(Apath, gt, ms, op, sc)
%DRIVER FOR SUBMODULAR OPTIMIZATION
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: November 16, 2023

disp('Job Starting');
disp(gt);
disp(ms);
disp(op);
disp(sc);
disp(Apath);

% Add paths
addpath(genpath('/home/jpic/DMD_gene/matlab2/'))

% Load data
load(Apath)

% Call optimizer
submodularSensorSelection(A,'gramT', gt, 'maxSensors', ms, 'outPath', op,'subCriteria', sc)

end