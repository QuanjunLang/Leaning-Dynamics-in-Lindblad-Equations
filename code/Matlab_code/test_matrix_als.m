clc
close all
clear all
rng(1)
% Matrix ALS test
% Suppose Bk = U*Ak*V;
% Estimate U and V from observations Ak and Bk

%% True para
n = 4;
m = 4;
K = 1000;

U_true = rand(m, m) + 1i*rand(m, m);
% U_true = U_true./U_true(1, 1);
V_true = U_true';



% V_true = rand(n, n);
A = rand(m, n, K);
temp = pagemtimes(U_true, A);
B = pagemtimes(temp, V_true);


%% Iteration
niter = 20;
all_U = zeros(m, m, niter);
all_V = zeros(n, n, niter);
all_V(:, :, 1) = randn(n, n);
V_cur = all_V(n, n, 1);

for i = 1:niter
    % Estimate U first 
    AA_U = cell(K, 1);
    bb = cell(K, 1);

    for k = 1:K
        AA_U{k} = kron((A(:, :, k)*V_cur).', eye(m, m));
        bb{k} = reshape(B(:, :, k), [], 1);
    end
    A_U = cat(1, AA_U{:});
    b_U = cat(1, bb{:});

    U_est = reshape(A_U\b_U, [m, m]);
    % U_est = U_est./U_est(1, 1);

    all_U(:, :, i) = U_est;
    U_cur = U_est;

    % % Estimate V next 
    AA_V = cell(K, 1);
    % bb = cell(K, 1);

    for k = 1:K
        AA_V{k} = kron(eye(n, n), U_cur*A(:, :, k));
        % bb{k} = reshape(B(:, :, k), [], 1);
    end

    A_V = cat(1, AA_V{:});
    b_V = cat(1, bb{:});

    V_est = reshape(A_V\b_V, [n, n]);

    % V_est = U_est;
    all_V(:, :, i) = V_est;
    V_cur = V_est;


    err_B(i) = norm(B - pagemtimes(pagemtimes(U_est, A), V_est), 'fro');
end


err_V = squeeze(sum((all_V - V_true).^2, [1,2]));
err_U = squeeze(sum((all_U - U_true).^2, [1,2]));


figure;hold on
plot(log10(err_V), ':', 'LineWidth', 3)
plot(log10(err_U))
plot(log10(err_B))

%%
U_est = all_U(:,  :, end);
V_est = all_V(:,  :, end);


B_est = pagemtimes(pagemtimes(U_est, A), V_est);

% err_B = 