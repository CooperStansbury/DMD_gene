function [Obft] = fbobsvt(Ab, Af, C, tb, tf)
%FBOBSVT obsv for a set number of time steps forward and backward
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 19, 2023

Of = obsvt(Af,C,tf);
Ob = flip(obsvt(Ab, flip(C,1), tb),1);
Of(1:size(C,1),:) = [];
Obft = [Ob; Of];

end

