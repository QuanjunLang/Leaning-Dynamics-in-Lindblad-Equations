function result = ALS_hamiltonian_kossakowski(all_rho_pair_data, sysInfo, trueInfo)

H = trueInfo.H_true;
v = trueInfo.v_true;
u = trueInfo.v_true;
F = trueInfo.F;
K = trueInfo.Kossakowski;
n = sysInfo.n;
p = sysInfo.p;

G = cell(n^2-1, n^2-1);
for k = 1:n^2-1
    for l = 1:n^2-1
        G_kl = kron(conj(F{k}), F{l}) - 0.5*kron((F{k}'*F{l}).', eye(n)) - 0.5*kron(eye(n), F{k}'*F{l});
        G{k, l} = sparse(G_kl);
    end
end


c = zeros(n^2-1, 1);
for j = 1:n^2-1
    c(j) = trace(H*F{j}');
end

tic
all_rho = squeeze(all_rho_pair_data(:, :, 1, :));
all_drho = squeeze(all_rho_pair_data(:, :, 2, :));
[~, n, M] = size(all_rho);
G_rho = cell(n^2-1, n^2-1, M);
for m = 1:M
    rho = all_rho(:, :, m);
    for k = 1:n^2-1
        for l = 1:n^2 - 1
            G_rho{k, l, m} = G{k, l}*vec(rho);
        end
    end
end
toc

%% DEBUG
% lambda = 1;
% tic;u_est = ALS_Kossakowski_u(all_rho_pair_data, v, H, G, p, lambda);toc
% tic;v_est = ALS_Kossakowski_v(all_rho_pair_data, u, H, G, p, lambda);toc
% tic;[H_est, c_est] = ALS_Kossakowski_H(all_rho_pair_data, u, v, G, F);toc
% K_est = u_est*v_est';



% sanity check
% norm(K - K_est)
% norm(u - u_est)
% norm(v - v_est)
% norm(c - c_est)

%% Iteration
niter = 400;

all_u = cell(niter, 1);
all_v = cell(niter, 1);
all_K = cell(niter, 1);
all_c = cell(niter, 1);
all_H = cell(niter, 1);

u0 = randn(size(u));
v0 = randn(size(v));

lambda = 0;
for i = 1:niter
    i
    % if i < 30iji
    %     lambda = 0;
    % else
    %     lambda = 0.01;
    % end

    [H0, c0] = ALS_Kossakowski_H(all_rho_pair_data, u0, v0, G, F);

    u0 = ALS_Kossakowski_u(all_rho_pair_data, v0, H0, G, p, lambda);

    v0 = ALS_Kossakowski_v(all_rho_pair_data, u0, H0, G, p, lambda);

    K0 = u0*v0';
    
    % K0 = (K0 + K0')/2;
    % [U, S, V] = svd(K0);
    % u0 = U(:, 1:p)*sqrt(S(1:p, 1:p));
    % v0 = V(:, 1:p)*sqrt(S(1:p, 1:p));

    all_u{i} = u0;
    all_v{i} = v0;
    all_K{i} = K0;
    all_c{i} = c0;
    all_H{i} = H0;
end



%%
u_err = zeros(niter, 1);
v_err = zeros(niter, 1);
K_err = zeros(niter, 1);
c_err = zeros(niter, 1);
H_err = zeros(niter, 1);

uv_df = zeros(niter, 1);
K_sym = zeros(niter, 1);

for i = 1:niter
    u_err(i) = norm(all_u{i} - u, 'fro');
    v_err(i) = norm(all_v{i} - v, 'fro');
    K_err(i) = norm(all_K{i} - K, 'fro');
    c_err(i) = norm(all_c{i} - c, 'fro');
    H_err(i) = norm(all_H{i} - H, 'fro');
    
    uv_df(i) = norm(all_u{i} - all_v{i}, 'fro');
    K_sym(i) = norm(all_K{i} - all_K{i}', 'fro');
end


lnwd = 1.5;
figure;hold on;
% plot(log10(u_err), '-o',  'DisplayName','u', 'LineWidth',lnwd);
% plot(log10(v_err), 'DisplayName','v', 'LineWidth',lnwd);
plot(log10(K_err), 'DisplayName','K', 'LineWidth',lnwd+3);
plot(log10(c_err), 'DisplayName','c', 'LineWidth',lnwd+3);
% plot(log10(H_err), 'DisplayName','H', 'LineWidth',lnwd);
% plot(log10(uv_df), 'DisplayName','u-v', 'LineWidth',lnwd);
% plot(log10(K_sym), 'DisplayName','K-K^t', 'LineWidth',lnwd);
legend()

result.u = all_u{end};
result.v = all_v{end};
result.K = all_K{end};
result.c = all_c{end};
result.H = all_H{end};
end