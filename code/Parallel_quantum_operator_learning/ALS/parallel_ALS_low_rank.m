function [X_est, outputInfo] = parallel_ALS_low_rank(A, b, r, varargin)
%% parallel_ALS_low_rank
% This function performs parallel reconstruction of low-rank matrices using
% Alternating Least Squares (ALS). The reconstruction is performed for sub-blocks,
% and the results are combined deterministically to form the final estimated matrix.
%
% Inputs:
%   - A: Tensor of size (n, n, d), representing the measurement operator
%   - b: Measurement vector
%   - r: Target rank for reconstruction
%   - varargin: Optional arguments (debugON, plotON, X_true, etc.)
%
% Outputs:
%   - X_est: Reconstructed large matrix
%   - outputInfo: Structure containing reconstruction details (errors, time, etc.)


% copyright - Quanjun Lang, 2024
%% Parse inputs using inputParser
p = inputParser;

% Required inputs
addRequired(p, 'A');
addRequired(p, 'b');
addRequired(p, 'r');

% Optional arguments with default values
addOptional(p, 'debugON', 0); 
addOptional(p, 'plotON', 0);

addOptional(p, 'X_true', 0);
addOptional(p, 'X_true_sub_blocks', 0);
addOptional(p, 'X_true_obs_blocks', 0);

addOptional(p, 'num_retry', 5);

% Parse arguments
parse(p, A, b, r, varargin{:});

% Extract parsed inputs
plotON              = p.Results.plotON;
debugON             = p.Results.debugON;
X_true_sub_blocks   = p.Results.X_true_sub_blocks;
X_true_obs_blocks   = p.Results.X_true_obs_blocks;
X_true              = p.Results.X_true;
num_retry           = p.Results.num_retry;

% Determine dimensions of the input
[N_o, M] = size(b); % Number of observed blocks
[n, ~, ~] = size(A); % Number of sub-blocks (assumes square structure)

% Initialize cell array for large matrix blocks
X_est_blocks = cell(n, n);
X_obs_blocks = cell(2*n-1, 1);
ALS_info = cell(n, 2);

% Start timer for performance evaluation
parallel_ALS_counter = tic;

% Start with extracting the measurements for the first row. 

new_b = zeros(2*n-1, M);

new_b(1:n, :) = b(1:n, :);

for i = 1:n-1
    new_b(i+n, :) = (b(n+i, :) -1i*b(2*n-1+i, :))/2;
end

%% Parallel reconstruction of the sub-blocks

parfor k = 1:2*n-1
    % k
    if k <= n
        % bk = b(k, :).';
        fprintf('ALS for the (%d, %d)-th\t block:', k, k);
    else
        % bk = b_fisrt_row(k-n, :).';
        fprintf('ALS for the (1, %d)-th\t block:', k-n+1);
    end
    
    for i = 1:num_retry
        [X_obs_blocks{k}, ALS_info{k}] = ALS(A, new_b(k, :).', r, 'NesterovON', 1);
        flag = ALS_info{k}.convergence_flag;
        if flag == 1
            break
        else
            fprintf('%d-th retry: ', i-1);
        end
    end
end



%% Deterministic reconstruction of the large matrix

% % Assign diagonal blocks from observed blocks
for i = 1:n
    X_est_blocks{i, i} = X_obs_blocks{i};
end
% 
% Assign off-diagonal blocks using complex-valued representation
for k = 2:n
    % Compute complex-valued blocks
    X_est_blocks{1, k} = X_obs_blocks{k+n-1};
    X_est_blocks{k, 1} = X_obs_blocks{k+n-1}';
end

% Reconstruct the remaining blocks using SVD and basis projections
for k = 2:n
    X_1k = X_est_blocks{1, k}; % Current off-diagonal block
    X_kk = X_est_blocks{k, k}; % Diagonal block for the k-th row/column
    
    % Perform SVD on the current off-diagonal block
    [U, S, V] = svd(X_1k);
    
    % Construct projection matrix C using diagonal block and basis
    C = X_kk * V(:, 1:r);
    
    % Compute remaining off-diagonal blocks
    for l = 2:n
        if l == k
            continue;
        end
        X_1l = X_est_blocks{1, l}; % Other off-diagonal block
        basis = S(1:r, 1:r) \ U(:, 1:r)' * X_1l; % Project onto SVD basis
        X_est_blocks{k, l} = C * basis; % Reconstruct block
    end
end


% Combine sub-blocks into the full reconstructed matrix
X_est = zeros(n^2, n^2);
for i = 1:n
    for j = 1:n
        X_est((i-1)*n+1:i*n, (j-1)*n+1:j*n) = X_est_blocks{i, j};
    end
end

% Store elapsed time
outputInfo.time = toc(parallel_ALS_counter);

%% Error analysis
% Compute errors for individual blocks
error_blocks = zeros(n, n);
for i = 1:n
    for j = 1:n
        error_blocks(i, j) = norm(X_est_blocks{i, j} - X_true_sub_blocks{i, j}, 'fro');
    end
end

% Compute Frobenius norms of true blocks
norm_blocks = zeros(n, n);
for i = 1:n
    for j = 1:n
        norm_blocks(i, j) = norm(X_true_sub_blocks{i, j}, 'fro');
    end
end

% Compute relative errors for individual blocks
rel_error_blocks = error_blocks ./ norm_blocks;

% Compute errors for the entire matrix
error_X = norm(X_est - X_true, 'fro');
rel_error_X = error_X / norm(X_true, 'fro');

% Store additional output information
outputInfo.ALS_info = ALS_info;
outputInfo.X_est_blocks = X_est_blocks;
outputInfo.error_blocks = error_blocks;
outputInfo.rel_error_blocks = rel_error_blocks;
outputInfo.norm_blocks = norm_blocks;
outputInfo.error_X = error_X;
outputInfo.rel_error_X = rel_error_X;
end
