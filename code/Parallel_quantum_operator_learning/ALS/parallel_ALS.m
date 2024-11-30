function [X_est, outputInfo] = parallel_ALS(A, b, r, varargin)






%%
p = inputParser;

addRequired(p, 'A');
addRequired(p, 'b');
addRequired(p, 'r');
addOptional(p, 'debugON', 0); 
addOptional(p, 'plotON', 0);
addOptional(p, 'X_true', 0);
addOptional(p, 'X_true_sub_blocks', 0);
addOptional(p, 'X_true_obs_blocks', 0);

parse(p, A, b, r, varargin{:});

plotON          = p.Results.plotON;
debugON         = p.Results.debugON;
X_true_sub_blocks    = p.Results.X_true_sub_blocks;
X_true_obs_blocks    = p.Results.X_true_obs_blocks;


[N_o, ~] = size(b);
[n, ~, ~] = size(A);

X_obs_blocks = cell(N_o, 1);
ALS_info = cell(N_o, 1);


%% Parallel reconstruction of the sub blocks
parfor k = 1:N_o
    if k <= n
        rk = r;
    else
        rk = 2*r;
    end
    bk = b(k, :).';
    fprintf('ALS for the %d-th \t block:', k)
    [X_obs_blocks{k}, ALS_info{k}] = ALS(A, bk, rk);
    % [X_obs_blocks{k}, ALS_info{k}] = ALS(A, bk, rk, 'maxIter', 800, 'X_true', X_true_obs_blocks{k}, 'debugON', 1, 'NesterovON', 0);
end

%% Deterministic reconstrurtion of the large matrix
X_est_blocks = cell(n, n);
for i = 1:n
    X_est_blocks{i, i} = X_obs_blocks{i};
end

for k = 2:n
    Rk = X_obs_blocks{k+n-1};
    Ik = X_obs_blocks{k+2*n-2};

    X_est_blocks{1, k} = (Rk - 1i*Ik)/2;
    X_est_blocks{k, 1} = (Rk + 1i*Ik)/2;
end


for k = 2:n
    X_1k = X_est_blocks{1, k};
    X_kk = X_est_blocks{k, k};
    [U, S, V] = svd(X_1k);
    C = X_kk * V(:, 1:r);
    for l = 2:n
        if l == k
            continue
        end
        X_1l = X_est_blocks{1, l};
        basis = S(1:r, 1:r)\U(:, 1:r)'*X_1l;
        X_est_blocks{k, l} = C * basis;
    end
end
    

X_est = zeros(n^2, n^2);
for i = 1:n
    for j = 1:n
        X_est((i-1)*n+1:i*n, (j-1)*n+1:j*n) = X_est_blocks{i, j};
    end
end

%% Error matrix
error_blocks = zeros(n, n);
for i = 1:n
    for j = 1:n
        error_blocks(i, j) = norm(X_est_blocks{i, j} - X_true_sub_blocks{i, j}, 'fro');
    end
end

norm_blocks = zeros(n, n);
for i = 1:n
    for j = 1:n
        norm_blocks(i, j) = norm(X_true_sub_blocks{i, j}, 'fro');
    end
end

rel_error_blocks = error_blocks ./ norm_blocks;



outputInfo.ALS_info = ALS_info;
outputInfo.X_est_blocks = X_est_blocks;
outputInfo.error_blocks = error_blocks;
outputInfo.rel_error_blocks = rel_error_blocks;
outputInfo.norm_blocks = norm_blocks;
end





% 
% %%
% X_22 = X_true_blocks{2, 2};
% X_12 = X_true_blocks{1, 2};
% X_23 = X_true_blocks{2, 3};
% X_13 = X_true_blocks{1, 3};
% 
% [U, S, V] = svd(X_12);
% 
% V(:, 1:r)' - S(1:r, 1:r)\U(:, 1:r)'*X_12        % The matrix {S(1:r, 1:r)\U(:, 1:r)'} project the X_12 block onto basis V(:, 1:r)'
% 
% 
% C = X_22 * V(:, 1:r);
% 
% C * V(:, 1:r)' - X_22                           % The matrix {C = X_22 * V(:, 1:r)} express X_22 by basis V(:, 1:r)'
% 
% % So the same setting works for other blocks
% 
% basis = S(1:r, 1:r)\U(:, 1:r)'*X_13;
% C * basis - X_23

