function [X_est, outputInfo] = parallel_ALS_first_row_new(A, b, r, varargin)
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
% p = inputParser;
% 
% % Required inputs
% addRequired(p, 'A');
% addRequired(p, 'b');
% addRequired(p, 'r');
% 
% % Optional arguments with default values
% addOptional(p, 'debugON', 0); 
% addOptional(p, 'plotON', 0);
% 
% addOptional(p, 'X_true', 0);
% addOptional(p, 'X_true_sub_blocks', 0);
% addOptional(p, 'X_true_obs_blocks', 0);
% 
% addOptional(p, 'num_retry', 0);
% 
% % Parse arguments
% parse(p, A, b, r, varargin{:});
% 
% % Extract parsed inputs
% plotON              = p.Results.plotON;
% debugON             = p.Results.debugON;
% X_true_sub_blocks   = p.Results.X_true_sub_blocks;
% X_true_obs_blocks   = p.Results.X_true_obs_blocks;
% X_true              = p.Results.X_true;
% num_retry           = p.Results.num_retry;
% 
% % Determine dimensions of the input
% [N_o, M] = size(b); % Number of observed blocks
% [n, ~, ~] = size(A); % Number of sub-blocks (assumes square structure)
% 
% % Initialize cell array for large matrix blocks
% X_est_blocks = cell(n, n);
% X_obs_blocks = cell(n-1, 1);
% ALS_info = cell(n, 2);
% 
% % Start timer for performance evaluation
% parallel_ALS_counter = tic;
% 
% % Start with extracting the measurements for the first row. 
% 
% new_b = zeros(n, M);
% 
% new_b(1, :) = b(1, :);
% 
% for i = 1:n-1
%     new_b(i+1, :) = (b(1+i, :) -1i*b(n+i, :))/2;
% end



%% Parse inputs using inputParser
p = inputParser;

% Required inputs
addRequired(p, 'A');
addRequired(p, 'b');
addRequired(p, 'r');

% Optional arguments with default values
addOptional(p, 'num_retry', 0);         % Number of retries when ALS does not converge


% addOptional(p, 'plotON', 0);
addOptional(p, 'debugON', 0);
addOptional(p, 'displayON', 0);
addOptional(p, 'compute_dnorm', 0);


addOptional(p, 'X_true', 0);            % True information to compare error
addOptional(p, 'X_true_sub_blocks', 0);

addOptional(p, 'Nesterov_beta', 0.5);
addOptional(p, 'ALS_rel_err_tol', 1e-8);
addOptional(p, 'ALS_loss_tol', 1e-8);
addOptional(p, 'ALS_debugON', 0)

% Parse arguments
parse(p, A, b, r, varargin{:});
% Extract parsed inputs
% plotON              = p.Results.plotON;
debugON             = p.Results.debugON;
displayON           = p.Results.displayON;
compute_dnorm       = p.Results.compute_dnorm;
X_true_sub_blocks   = p.Results.X_true_sub_blocks;
X_true              = p.Results.X_true;
num_retry           = p.Results.num_retry;

Nesterov_beta        = p.Results.Nesterov_beta;
ALS_rel_err_tol     = p.Results.ALS_rel_err_tol;
ALS_loss_tol        = p.Results.ALS_loss_tol;
ALS_debugON         = p.Results.ALS_debugON;
%% Extract parameters from input matrices

fprintf('\nLearning operator in parallel for all first row blocks:\n')

% Determine dimensions of the input
[~, M] = size(b); % Number of observed blocks
[N, ~, ~] = size(A); % Number of sub-blocks (assumes square structure)

% Start timer for performance evaluation
time_counter = tic;

% Start with extracting the measurements for the first row.
new_b = zeros(N, M);
new_b(1, :) = b(1, :);
for i = 1:N-1
    new_b(i+1, :) = (b(1+i, :) -1i*b(N+i, :))/2;
end

% new_b serves are the measurement result as if we take observables to be
% E_11, ..., E_1N.



% Initialize cell array for large matrix blocks
X_est_blocks = cell(N, N);
X_obs_blocks = cell(N-1, 1);
ALS_info = cell(N, 2);


%% Parallel reconstruction of the sub-blocks

for k = 1:N
    fprintf('ALS for the (1, %d)-th\t block:', k);
    [X_obs_blocks{k}, ALS_info{k}] = ALS(A, new_b(k, :).', r, 'NesterovON', 1, 'X_true', X_true_sub_blocks{1, k}, ...
        'Nesterov_beta', Nesterov_beta, 'loss_tol', ALS_loss_tol, 'rel_err_tol', ALS_rel_err_tol, 'debugON', ALS_debugON);
    flag = ALS_info{k}.convergence_flag;
    if flag ~= 1 && num_retry>0
        for i = 1:num_retry
            fprintf('%d-th retry: ', i);
            [X_obs_blocks{k}, ALS_info{k}] = ALS(A, new_b(k, :).', r, 'NesterovON', 1, 'X_true', X_true_sub_blocks{1, k}, ...
                'Nesterov_beta', Nesterov_beta, 'loss_tol', ALS_loss_tol, 'rel_err_tol', ALS_rel_err_tol, 'debugON', ALS_debugON);
            flag = ALS_info{k}.convergence_flag;
            if flag == 1
                break
            end
        end
    end
    if flag ~= 1
        fprintf('Recovery for this block failed. Please consider generate more data, or increase the retry numbers\n');
        X_est = [];
        outputInfo.rel_error_fro_X = 100;
        outputInfo.flag = flag;
        outputInfo.time = ALS_info{k}.time;
        outputInfo.ALS_info = ALS_info{k};
        return 
    end
end

X_est_first_row_blocks = zeros(N, N^2);
for i = 1:N
    X_est_first_row_blocks(:, (i-1)*N+1:i*N) = X_obs_blocks{i};
end


%% Deterministic reconstruction of the large matrix

% % Assign diagonal blocks from observed blocks
for k = 1:N
    X_est_blocks{1, k} = X_obs_blocks{k};
    X_est_blocks{k, 1} = X_obs_blocks{k}';
end

X_11 = X_est_blocks{1, 1}; % Current off-diagonal block
% Perform SVD on the current off-diagonal block
[U, S, V] = svd(X_11);

% Reconstruct the remaining blocks using SVD and basis projections
for k = 2:N
    
    X_k1 = X_est_blocks{k, 1}; % Diagonal block for the k-th row/column
   
    % Construct projection matrix C using diagonal block and basis
    C = X_k1 * V(:, 1:r);
    
    % Compute remaining off-diagonal blocks
    for l = 2:N
        X_1l = X_est_blocks{1, l}; % Other off-diagonal block
        basis = S(1:r, 1:r) \ U(:, 1:r)' * X_1l; % Project onto SVD basis
        X_est_blocks{k, l} = C * basis; % Reconstruct block
    end
end


% Combine sub-blocks into the full reconstructed matrix
X_est = zeros(N^2, N^2);
for i = 1:N
    for j = 1:N
        X_est((i-1)*N+1:i*N, (j-1)*N+1:j*N) = X_est_blocks{i, j};
    end
end


%% Store estimation results
outputInfo.time = toc(time_counter);
outputInfo.ALS_info = ALS_info;
outputInfo.flag = flag;


outputInfo.X_est_blocks = X_est_blocks;
% outputInfo.X_est_first_row_blocks = X_est_first_row;

fprintf('Finished in %.2f seconds \n', outputInfo.time)
%% Error analysis


outputInfo = error_analysis(outputInfo, X_est, X_true, 'X_est_blocks', X_est_blocks, 'X_true_sub_blocks', X_true_sub_blocks, 'displayON', displayON, 'compute_dnorm', compute_dnorm);



% if debugON && ismatrix(X_true)
% 
%     % Frobenius norm for each submatrix
%     error_fro_blocks = zeros(N, N);
%     for i = 1:N
%         for j = 1:N
%             error_fro_blocks(i, j) = norm(X_est_blocks{i, j} - X_true_sub_blocks{i, j}, 'fro');
%         end
%     end
% 
%     % Frobenius norms of true blocks
%     norm_blocks = zeros(N, N);
%     for i = 1:N
%         for j = 1:N
%             norm_blocks(i, j) = norm(X_true_sub_blocks{i, j}, 'fro');
%         end
%     end
% 
%     % Relative Frobenius norm of each submatrix
%     rel_error_fro_blocks = error_fro_blocks ./ norm_blocks;
% 
%     % Frobenius errors for the entire matrix
%     error_fro_X = norm(X_est - X_true, 'fro');
%     rel_error_fro_X = error_fro_X / norm(X_true, 'fro');
% 
% 
%     outputInfo.error_fro_blocks = error_fro_blocks;
%     outputInfo.rel_error_fro_blocks = rel_error_fro_blocks;
%     outputInfo.norm_blocks = norm_blocks;
%     outputInfo.error_fro_X = error_fro_X;
%     outputInfo.rel_error_fro_X = rel_error_fro_X;
% 
% 
%     % Nuclear norm error of reshaped matrix X
%     error_nuc_X     = sum(svd(X_est - X_true));
%     rel_error_nuc_X = error_nuc_X./sum(svd(X_true));
% 
%     outputInfo.error_nuc_X = error_nuc_X;
%     outputInfo.rel_error_nuc_X = rel_error_nuc_X;
% 
% 
%     % Nuclear norm error of operator
%     error_nuc_op     = sum(svd(rearrangement_R(X_est) - rearrangement_R(X_true)));
%     rel_error_nuc_op = error_nuc_op./sum(svd(rearrangement_R(X_true)));
% 
%     outputInfo.error_nuc_op = error_nuc_op;
%     outputInfo.rel_error_nuc_op = rel_error_nuc_op;
% 
% 
% end


%% Ploting and printing
if outputInfo.flag && debugON && ismatrix(X_true)
    figure;
    subplot(131)
    % surf(RL_Info.error_blocks);view(0, 90);
    imagesc(outputInfo.error_fro_blocks)
    colorbar
    title('L blocks error')
    subplot(132)
    % surf(RL_Info.rel_error_blocks);view(0, 90);
    imagesc(outputInfo.rel_error_fro_blocks)
    colorbar
    title('L blocks relative error')
    subplot(133)
    % surf(RL_Info.norm_blocks);view(0, 90);
    imagesc(outputInfo.norm_blocks)
    colorbar
    title('L blocks norm')
end


% if debugON
%     fprintf('Relative Frobenius error \t= %.8f, \tRelative Nuclear error on Choi \t= %.8f, \tRelative Nuclear error on operator \t= %.8f \n',rel_error_fro_X, rel_error_nuc_X, rel_error_nuc_op);
%     fprintf('Frobenius error \t\t= %.8f, \tNuclear error on Choi \t\t= %.8f, \tNuclear error on operator \t\t= %.8f \n',error_fro_X, error_nuc_X, error_nuc_op);
% end




% 
% % Store elapsed time
% outputInfo.time = toc(time_counter);
% 
% %% Error analysis
% % Compute errors for individual blocks
% error_blocks = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         error_blocks(i, j) = norm(X_est_blocks{i, j} - X_true_sub_blocks{i, j}, 'fro');
%     end
% end
% 
% % Compute Frobenius norms of true blocks
% norm_blocks = zeros(N, N);
% for i = 1:N
%     for j = 1:N
%         norm_blocks(i, j) = norm(X_true_sub_blocks{i, j}, 'fro');
%     end
% end
% 
% % Compute relative errors for individual blocks
% rel_error_blocks = error_blocks ./ norm_blocks;
% 
% % Compute errors for the entire matrix
% error_X = norm(X_est - X_true, 'fro');
% rel_error_X = error_X / norm(X_true, 'fro');
% 
% % Store additional output information
% outputInfo.ALS_info = ALS_info;
% outputInfo.X_est_blocks = X_est_blocks;
% outputInfo.error_blocks = error_blocks;
% outputInfo.rel_error_blocks = rel_error_blocks;
% outputInfo.norm_blocks = norm_blocks;
% outputInfo.error_X = error_X;
% outputInfo.rel_error_X = rel_error_X;
% outputInfo.flag = flag;
% outputInfo.X_est_first_row_blocks = X_est_first_row_blocks;
% 
% 
% % Nuclear norm error 
% error_X_nuclear     = sum(svd(X_est - X_true));
% rel_error_X_nuclear = error_X_nuclear./sum(svd(X_true));
% 
% outputInfo.error_X_nuclear = error_X_nuclear;
% outputInfo.rel_error_X_nuclear = rel_error_X_nuclear;
% 
% 
% % Nuclear norm error of Choi
% error_X_nuclear     = sum(svd(X_est - X_true));
% rel_error_X_nuclear = error_X_nuclear./sum(svd(X_true));
% 
% outputInfo.error_X_nuclear = error_X_nuclear;
% outputInfo.rel_error_X_nuclear = rel_error_X_nuclear;
% 
% 
% % Nuclear norm error of operator
% error_X_nuclear_op     = sum(svd(rearrangement_R(X_est) - rearrangement_R(X_true)));
% rel_error_X_nuclear_op = error_X_nuclear_op./sum(svd(rearrangement_R(X_true)));
% 
% outputInfo.error_X_nuclear_op = error_X_nuclear_op;
% outputInfo.rel_error_X_nuclear_op = rel_error_X_nuclear_op;
% 
% end
