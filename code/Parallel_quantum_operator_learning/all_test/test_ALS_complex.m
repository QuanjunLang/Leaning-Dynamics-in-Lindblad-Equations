% clc
% close all
% clear all
% rng(10);
% 
% r = 3;
% N = 10;
% M = 80;
% sigma = 0.0001;
% 
% U = randn(N, r) + 1i*randn(N, r);
% V = randn(N, r) + 1i*randn(N, r);
% X = U*V';
% 
% A = zeros(N, N, M);
% b = zeros(N, 1);
% for m = 1:M
%     A(:, :, m) = randn(N, N) + 1i*randn(N, N);
%     b(m) = sum(conj(A(:, :, m)).*X, 'all') + randn*sigma;
% end
% 
% [M_est, outputInfo] = ALS(A, b, r)
% all_error = sqrt(squeeze(sum(abs(outputInfo.all_X - X).^2, [1,2])));
% 
% figure;
% plot(log10(all_error))
% 
