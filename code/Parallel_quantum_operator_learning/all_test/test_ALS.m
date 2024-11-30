clc
close all
clear all
rng(1)
addPaths
%%
M = 50;
N = 5;
r = 3;


%%
% Generate random sensing matrices A_i
A = cell(M, 1); % Cell array to store sensing matrices
for m = 1:M
    A{m} = randn(N, N) + randn(N, N)*1i; % Random sensing matrix
end

% Generate ground truth matrix M_true of rank r
U_true = randn(N, r) + randn(N, r)*1i;
V_true = randn(N, r) + randn(N, r)*1i;
M_true = U_true * V_true';

% Generate observations b_i = <A_i, M_true>
b = zeros(M, 1);
for m = 1:M
    b(m) = sum(conj(A{m}) .* M_true, 'all'); % Frobenius inner product
end



%%
M_est = ALS(A, b, r);
M_true