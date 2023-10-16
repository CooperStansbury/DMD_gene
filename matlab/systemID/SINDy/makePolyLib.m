function lib = makePolyLib(n, order)
%MAKEPOLYLIB Makes a library of polynomial functions
%
%   INPUTS:
%       - n: number of variables
%       - order: polynomial order
%   OUTPUTS:
%       - lib: cell array of polynomial functions
%
%{
    D = rand(3,10);
    F = firstDifferences(D,1); F = F{1};
    lib = { @(D)D(1,:),
            @(D)D(2,:),
            @(D)D(3,:),
            @(D)D(1,:).*D(2,:),
            @(D)D(1,:).*D(3,:),
            @(D)D(2,:).*D(3,:)  };
%}
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 2, 2023


if order ~= 3
    error("DONT DO THIS");
end

lib = cell(nchoosek(n,order), 1);
vars = nchoosek(1:n,order);
for i=1:length(lib)
    f = "@(D)D(" + string(vars(i,1)) + ",:).*D(" + string(vars(i,2)) + ",:).*D(" + string(vars(i,3)) + ",:);";
    f = eval(f);
    lib{i} = f;
    % @(D)D(vars(i,1),:).*D(vars(i,2),:).*D(vars(i,3),:);
end

end