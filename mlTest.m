load('HW6-1b.mat'); 
% ML 最大似然检测
% Map x_hat_noisy to nearest 16-QAM constellation
qam_constellation = [-3-3j, -3-1j, -3+1j, -3+3j, ...
                     -1-3j, -1-1j, -1+1j, -1+3j, ...
                      1-3j,  1-1j,  1+1j,  1+3j, ...
                      3-3j,  3-1j,  3+1j,  3+3j];

[X1, X2, X3] = ndgrid(qam_constellation, qam_constellation, qam_constellation);
X_all = [X1(:), X2(:), X3(:)].';

num_candidates = size(X_all, 2);
Gamma = zeros(1, num_candidates);  

for i = 1:num_candidates
    x_candidate = X_all(:, i); 
    Gamma(i) = norm(yprime - Hmatrix * x_candidate)^2;
end

[Gamma_min, min_idx] = min(Gamma);
x_ML = X_all(:, min_idx);

figure;
plot(1:num_candidates, Gamma, '-o');
xlabel('Index of Candidate');
ylabel('Gamma(x)');
title('Cost Function vs Candidate Index');
grid on;

disp(['Minimum cost function value: ', num2str(Gamma_min)]);
disp('ML solution: ');
disp(x_ML);
