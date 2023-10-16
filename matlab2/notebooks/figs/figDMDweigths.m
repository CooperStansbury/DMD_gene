function figDMDweights(A)

AA = abs(A(:));
SA = sort(AA,'descend');

k = 1e4;

figure;
scatter(1:k,SA(1:k),'.');
set(gca, 'YScale', 'log');
ylabel('Absolute Value');
xlabel('Top Interactions');
title('DMD Weights');

end
