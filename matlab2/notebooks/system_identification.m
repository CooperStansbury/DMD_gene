%% System Identification
%
%   This notebook explores alternative methods for approximating cell
%   dynamics. DMD, SINDy, and other data driven system identification
%   frameworks often exploit a form of the regression:
%
%                          min_C ||A-B*C||
%
%   where A and B are matrices derived from time series data. However, 1) 
%   for high dimensional or nonlinear problems the regression can become
%   computationally intractible, and 2) alternative data modalities beyond
%   time series data may be helpful for the system identification problem.
%   We can use structural, static, or other data modalities to identify the
%   sparsity of C prior to performing the regression. This is formulated as
%   the constrained optimization:
%
%               min_C ||A-B*C|| such that (C ~= 0) == D
%
%   where D is derived from alternative data collected using orthognal
%   biological assays.
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: September 18, 2023

%% structuredLS
%   The structuredLS function solves the minimization problem:
%
%               min_C ||A-B*C|| such that (C ~= 0) == D
%
%   To demonstrate this, matrices A, B, Ctrue, and D are constructed. Then,
%   A,B, and D are passed as arguments to structuredLS, which produces
%   Clearned, and Cdense is obtained as the unconstrained least squares
%   solution. Clearned and Cdense are compared in terms of their
%   approximation to Ctrue and the time it takes to compute each. We see
%   that Clearned has a better approximation error and Cdense is faster to
%   compute.

clear

n = 300;
m = 400;
p = m;

B = rand(n,p);
Ctrue = full(sprand(p,m,0.2));
S = (Ctrue ~= 0);
A = B * Ctrue;

tic;
Clearned = structuredLS(A,B,S);
ts = toc;
tic;
Cdense = B\A;
td = toc;

nd = norm(Ctrue-Cdense);
ns = norm(Ctrue-Clearned);

disp("A: " + string(size(A)));
disp("B: " + string(size(B)));
disp("C: " + string(size(Ctrue)));
disp("Sparse Error: " + string(ns) + ", Time: " + string(ts));
disp("Dense Error: " + string(nd) + ", Time: " + string(td));


%% Can structuredLS compete with DMD?
%
%   In this section, I take the above results and attempt to use the same
%   structuredLS problem for approximating sparse linear systems.
%   Currently, classic DMD performs better than structuredLS, but I have
%   not yet figured out why.

clear;

n = 200;
t = 1000;

Atrue = full(sprand(n,n,0.05));
Atrue = Atrue ./ sum(Atrue,1);
% sum(Atrue,1)
X = zeros(n,t);
X(:,1) = rand(n,1);
for i=2:t
    X(:,i) = Atrue * X(:,i-1);
end

disp(find(isnan(X)));


% t = 20000;
% X = X(:,1:t);

tic;
out = DMD(X,[],1);
td = toc; disp('DMD');
tic;
As = SDMD(X,(Atrue ~= 0));
ts = toc; disp('Sparse');
tic;
Ae = EDMD(X);
te = toc; disp('Exact');

nd = norm(Atrue - out.DMD.A_bar);
ns = norm(Atrue - As);
ne = norm(Atrue - Ae);

X1 = X(:,1:end-1);
X2 = X(:,2:end);

Xd2 = out.DMD.A_bar * X1;
ed = sum(abs(X2 - Xd2),'all');
Xs2 = As * X1;
es = sum(abs(X2 - Xs2),'all');
Xe2 = Ae * X1;
ee = sum(abs(X2 - Xe2),'all');

disp("X: " + string(size(X)));
disp("A: " + string(size(Atrue)));
disp("Sparse Error: " + string(ns) + ", Time: " + string(ts) + ", Prediction Error: " + string(es));
disp("DMD Error: " + string(nd) + ", Time: " + string(td) + ", Prediction Error: " + string(ed));
disp("Exact Error: " + string(ne) + ", Time: " + string(te) + ", Prediction Error: " + string(ee));

