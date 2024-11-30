function result = L_decomposition_hamiltonian_kossakowski(L, p, trueInfo, varargin)

parser = inputParser;

addRequired(parser, 'L');
addRequired(parser, 'p');
addRequired(parser, 'trueInfo');
addOptional(parser, 'niter', 50);
addOptional(parser, 'threshold', 1e-8);
addOptional(parser, 'ttl', '');
addOptional(parser, 'plotON', 1);
parse(parser, L, p, trueInfo, varargin{:});

niter           = parser.Results.niter;
L               = parser.Results.L;
p               = parser.Results.p;
trueInfo        = parser.Results.trueInfo;
threshold       = parser.Results.threshold;
ttl             = parser.Results.ttl;
plotON          = parser.Results.plotON;
%%
H = trueInfo.H_true;
L_vec = vec(L);
v = trueInfo.v_true;
u = trueInfo.v_true;
F = trueInfo.F;
K = trueInfo.Kossakowski;
n = sqrt(length(L));
h = trueInfo.h;
G = trueInfo.G;

% G = cell(n^2-1, n^2-1);
% for k = 1:n^2-1
%     for l = 1:n^2-1
%         G_kl = kron(conj(F{k}), F{l}) - 0.5*kron((F{k}'*F{l}).', eye(n)) - 0.5*kron(eye(n), F{k}'*F{l});
%         G{k, l} = sparse(G_kl);
%     end
% end
%
%
% c = zeros(n^2-1, 1);
% for j = 1:n^2-1
%     c(j) = trace(H*F{j}');
% end




%% DEBUG
% lambda = 0.00000001;
%
% tic;u_est = L_decomposition_Kossakowski_u(L_vec, v, H, G, n, p, lambda);toc
% tic;v_est = L_decomposition_Kossakowski_v(L_vec, u, H, G, n, p, lambda);toc
% tic;[H_est, c_est] = L_decomposition_Kossakowski_H(L, u, v, G, F, n);toc
% K_est = u_est*v_est';



% sanity check
% norm(K - K_est)
% norm(u - u_est)
% norm(v - v_est)
% norm(c - c_est)

%% Iteration
tic;

all_u = cell(2, 1);
all_v = cell(2, 1);
all_K = cell(2, 1);
all_h = cell(2, 1);
all_H = cell(2, 1);

u0 = randn(size(u));
v0 = randn(size(v));

lambda = 0;
for i = 1:niter
    % if i < 30
    %     lambda = 0;
    % else
    %     lambda = 0.01;
    % end
    % i
    [H0, h0] = L_decomposition_Kossakowski_H(L, u0, v0, G, F, n);
    u0 = L_decomposition_Kossakowski_u(L_vec, v0, H0, G, n, p, lambda);
    v0 = L_decomposition_Kossakowski_v(L_vec, u0, H0, G, n, p, lambda);
    K0 = u0*v0';

    % K0 = (K0 + K0')/2;
    % [U, S, V] = svd(K0);
    % u0 = U(:, 1:p)*sqrt(S(1:p, 1:p));
    % v0 = V(:, 1:p)*sqrt(S(1:p, 1:p));

    all_u{i} = u0;
    all_v{i} = v0;
    all_K{i} = K0;
    all_h{i} = h0;
    all_H{i} = H0;

    if i > 1
        c_diff = norm(h0 - all_h{i-1});
        K_diff = norm(K0 - all_K{i-1});
        if c_diff<threshold && K_diff < threshold
            niter = i;
            break
        end
    end
end
total_time = toc;


%%
u_err = zeros(2, 1);
v_err = zeros(2, 1);
K_err = zeros(2, 1);
h_err = zeros(2, 1);
H_err = zeros(2, 1);

uv_df = zeros(2, 1);
K_sym = zeros(2, 1);

for i = 1:niter
    u_err(i) = norm(all_u{i} - u, 'fro');
    v_err(i) = norm(all_v{i} - v, 'fro');
    K_err(i) = norm(all_K{i} - K, 'fro');
    h_err(i) = norm(all_h{i} - h, 'fro');
    H_err(i) = norm(all_H{i} - H, 'fro');

    uv_df(i) = norm(all_u{i} - all_v{i}, 'fro');
    K_sym(i) = norm(all_K{i} - all_K{i}', 'fro');
end

if plotON
    lnwd = 3;
    figure;hold on;
    % plot(log10(u_err), '-o',  'DisplayName','u', 'LineWidth',lnwd);
    % plot(log10(v_err), 'DisplayName','v', 'LineWidth',lnwd);
    plot(log10(K_err), 'DisplayName','K', 'LineWidth',lnwd);
    plot(log10(h_err), 'DisplayName','c', 'LineWidth',lnwd);
    % plot(log10(H_err), 'DisplayName','H', 'LineWidth',lnwd);
    % plot(log10(uv_df), 'DisplayName','u-v', 'LineWidth',lnwd);
    % plot(log10(K_sym), 'DisplayName','K-K^t', 'LineWidth',lnwd);
    legend()
    title(ttl);

end
result.u = all_u{end};
result.v = all_v{end};
result.K = all_K{end};
result.c = all_h{end};
result.H = all_H{end};

result.h_err = h_err;
result.K_err = K_err;

result.C = cell(p, 1);
[U, S, ~] = svd(result.K);
US = U*sqrt(S);
for i = 1:p
    result.C{i} = zeros(n, n);
    for j = 1:n^2-1
        result.C{i} = result.C{i} + US(j, i)*F{j};
    end
end

result.time = total_time;
end