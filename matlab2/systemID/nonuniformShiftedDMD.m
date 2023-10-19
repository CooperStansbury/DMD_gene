function [Output] = nonuniformShiftedDMD(D,MetaData,thresh)
%SHIFTEDDMD Summary of this function goes here
%   Detailed explanation goes here
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 18, 2023
%
% Code based on Joshua Proctor's at the Institute for Disease Modeling 2015
% with modifications to build the initial and final data matrices with
% multiple replicates and shifting according the Hasnain et. al. 2023

% Construct data matrices note, time snapshots are columns
replicates = numel(D);
X = [];
Xp = [];
for r=1:replicates
    DM = D{r}';
    X = [X DM(:,1:end-1)];
    Xp = [Xp DM(:,2:end)];
end

% Compute the SVD 
    [UX,SigX,VX] = svd(X,'econ');
    
% Choose the number of singular values based on energy percent
    r = find(cumsum((diag(SigX)/sum(diag(SigX)))) > thresh,1);
    
% Compute A_tilde    
    A = UX(:,1:r)'*Xp*VX(:,1:r)*inv(SigX(1:r,1:r));
    
% Compute Eigenanalysis
    [V,D] = eig(A);

% Compute Dynamic Modes
    Dyn = Xp*VX(:,1:r)*inv(SigX(1:r,1:r))*V;

% Place computed variables in a single output file 
Output.DataMatrix = D;
Output.MetaData   = MetaData;

Output.X = X;
Output.Xp   = Xp;
     
Output.DMD.D = diag(D);
Output.DMD.V = V;
Output.DMD.DynamicModes = Dyn;
Output.DMD.A = A;
Output.DMD.Sig = SigX;
Output.DMD.r = r;
Output.DMD.UX = UX;
Output.DMD.VX = VX;         

% if size(DataMatrix,1) <= 1000
% Compute A_bar SCOTT ADD
    A_bar = Xp*VX*pinv(SigX)*UX';
    [V_bar,D_bar] = eig(A_bar);
    
    Output.DMD.A_bar = A_bar;           %SCOTT ADD
    Output.DMD.D_bar = diag(D_bar);     %SCOTT ADD
    Output.DMD.V_bar = V_bar;           %SCOTT ADD
% end


end

