function [Ab, Af] = fbDMD(D,reps,t,thresh)
%FBDMD Forward-Backward Dynamic Mode Decomposition
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 18, 2023

% clear; % D = rand(4,12); % reps = 3; % T = 4; % t = 2;

T = size(D,2) / reps;

% forward data
Df = [];
for r=1:reps
    Df = [Df D(:,(r-1)*T+t+1:r*T)];
end

% backward data
Db = [];
for r=1:reps
    Db = [Db D(:,(r-1)*T+1:(r-1)*T+t)];
end
Db = flip(Db,2);

% shifted dmd
outf = shiftedDMD(Df, reps, [], thresh);
outb = shiftedDMD(Db, reps, [], thresh);
Af = outf.DMD.A_bar;
Ab = outb.DMD.A_bar;

end

