function pred = predictLinear(slm, signalsTest)
%PREDICTLINEAR Predict velocity with simple linear model
%
% Auth: Joshua Pickard
%       jpic@umich.edu
% Date: October 7, 2023

pred = cell(numel(signalsTest), 1);
for i=1:numel(pred)
    signals = signalsTest{i};
    pred{i} = slm' * signals';
end

end
