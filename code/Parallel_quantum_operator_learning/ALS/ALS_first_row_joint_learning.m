function [X_est, outputInfo] = ALS_first_row_joint_learning(A, b, r, varargin)
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
% %% Parse inputs using inputParser
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
% addOptional(p, 'num_retry', 0);
% 
% % Parse arguments
% parse(p, A, b, r, varargin{:});
% 
% % Extract parsed inputs
% plotON              = p.Results.plotON;
% debugON             = p.Results.debugON;
% X_true_sub_blocks   = p.Results.X_true_sub_blocks;
% X_true              = p.Results.X_true;
% num_retry           = p.Results.num_retry;
% 
% % Determine dimensions of the input
% [~, M] = size(b); % Number of observed blocks
% [n, ~, ~] = size(A); % Number of sub-blocks (assumes square structure)
% 
% % Initialize cell array for large matrix blocks
% X_est_blocks = cell(n, n);
% 
% % Start timer for performance evaluation
% time_counter = tic;
% 
% % Start with extracting the measurements for the first row. 
% new_b = zeros(n, M);
% new_b(1, :) = b(1, :);
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


addOptional(p, 'plotON', 0);
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
plotON              = p.Results.plotON;
debugON             = p.Results.debugON;
displayON           = p.Results.displayON;
compute_dnorm       = p.Results.compute_dnorm;
X_true_sub_blocks   = p.Results.X_true_sub_blocks;
X_true              = p.Results.X_true;
num_retry           = p.Results.num_retry;

Nesterov_beta       = p.Results.Nesterov_beta;
ALS_rel_err_tol     = p.Results.ALS_rel_err_tol;
ALS_loss_tol        = p.Results.ALS_loss_tol;
ALS_debugON         = p.Results.ALS_debugON;
%% Extract parameters from input matrices

fprintf('\nLearning operator from the first row jointly:\n')

% Determine dimensions of the input
[~, M] = size(b); % Number of observed blocks
[N, ~, ~] = size(A); % Number of sub-blocks (assumes square structure)

% Initialize cell array for large matrix blocks
X_est_blocks = cell(N, N);

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



%% Reconstruction of the first row
fprintf('ALS for the first row : ')

% ALS_debugON = 0;
NesterovON =1;

[X_est_first_row, ALS_outputInfo] = ALS_shared_U(A, new_b, r, 'NesterovON', NesterovON, 'debugON', ALS_debugON, 'X_true', X_true(1:N, :), 'Nesterov_beta', Nesterov_beta, ...
    'loss_tol', ALS_loss_tol, 'rel_err_tol', ALS_rel_err_tol);

flag = ALS_outputInfo.convergence_flag;
if flag ~= 1 && num_retry>0
    for i = 1:num_retry
        fprintf('%d-th retry: ', i);
        [X_est_first_row, ALS_outputInfo] = ALS_shared_U(A, new_b, r, 'NesterovON', NesterovON, 'debugON', ALS_debugON, 'X_true', X_true(1:N, :), 'Nesterov_beta', Nesterov_beta, ...
            'loss_tol', ALS_loss_tol, 'rel_err_tol', ALS_rel_err_tol);
        flag = ALS_outputInfo.convergence_flag;
        if flag == 1
            break
        end
    end
end
if flag ~= 1
    fprintf('Recovery for the first row failed. Please consider generate more data, or increase the retry numbers\n');
    X_est = [];
    outputInfo.rel_error_X = 100;
    outputInfo.flag = flag;
        outputInfo.time = ALS_outputInfo.time;
        outputInfo.ALS_info = ALS_outputInfo;
    return 
end

%% Deterministic reconstruction of the large matrix
X_est_rows = cell(N, 1);
X_est_rows{1} = X_est_first_row;

[~, ~, V] = svd(X_est_first_row);
V0 = V(:, 1:r);
V0_11 = V(1:N, 1:r);

for k = 2:N
    X_k1 = X_est_first_row(:, (k-1)*N+1:k*N)'; % Diagonal block for the k-th row/column
    C = (pinv(V0_11)*(X_k1'))';
    X_est_rows{k} = C * V0';
end

% Reconstruct the remaining blocks using SVD and basis projections
for k = 1:N
    for l = 1:N
        X_est_blocks{k, l} = X_est_rows{k}(:, (l-1)*N+1:l*N);
    end
end


% Combine sub-blocks into the full reconstructed matrix
X_est = zeros(N^2, N^2);
for i = 1:N
    for j = 1:N
        X_est((i-1)*N+1:i*N, (j-1)*N+1:j*N) = X_est_blocks{i, j};
    end
end

% % Store elapsed time
% outputInfo.time = toc(parallel_ALS_counter);
% 





%% Store estimation results
outputInfo.time = toc(time_counter);
outputInfo.ALS_info = ALS_outputInfo;
outputInfo.flag = flag;


outputInfo.X_est_blocks = X_est_blocks;
outputInfo.X_est_first_row_blocks = X_est_first_row;

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
if outputInfo.flag && debugON && ismatrix(X_true) && plotON
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

end




% %% Error analysis
% % Compute errors for individual blocks
% error_blocks = zeros(n, n);
% for i = 1:n
%     for j = 1:n
%         error_blocks(i, j) = norm(X_est_blocks{i, j} - X_true_sub_blocks{i, j}, 'fro');
%     end
% end
% 
% % Compute Frobenius norms of true blocks
% norm_blocks = zeros(n, n);
% for i = 1:n
%     for j = 1:n
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
% 
% 
% % Store additional output information
% outputInfo.ALS_info = outputInfo;
% outputInfo.X_est_blocks = X_est_blocks;
% outputInfo.error_blocks = error_blocks;
% outputInfo.rel_error_blocks = rel_error_blocks;
% outputInfo.norm_blocks = norm_blocks;
% outputInfo.error_X = error_X;
% outputInfo.rel_error_X = rel_error_X;
% outputInfo.flag = flag;
% outputInfo.X_est_first_row_blocks = X_est_first_row;
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
% end
