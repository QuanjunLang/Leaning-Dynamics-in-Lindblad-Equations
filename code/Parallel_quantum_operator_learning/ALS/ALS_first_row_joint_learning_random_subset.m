function [X_est, outputInfo] = ALS_first_row_joint_learning_random_subset(A, b, r, varargin)
%%
% This function learns the reshaped quantum operator from observation of
% sub-blocks in its first row.
%
% Suppose X is a quantum operator, reshaped to a matrix with (N*N)x(N*N).
% It is the Choi matrix of the channel operator in the quantum process
% tomography. We decompose X to N x N submatrices of size N*N.
%
% We assume that X has rank r << N^2. We also assume that we are able to
% observe the submatrices on the first row, denoted by X_11, ..., X_1N.
%
% By the low rank assumption, each of the submatrices will have rank at
% most r. We have X_1k = U*V_k', where U is shared accross each matrix.
%
%%%%%%%%%%%%%%% Learning Process: %%%%%%%%%%%%%
% 1. Random Selection and ALS Reconstruction:
%    - Randomly select m = ln(N) indices k1, ..., km from {1, ..., N}.
%    - Jointly reconstruct the selected blocks X_1k1, ..., X_1km using
%      Alternating Least Squares (ALS), assuming a shared U in the
%      decomposition X_1k = U * V_k'.
%
% 2. Least Squares Estimation:
%    - Use the estimated U from ALS to solve for remaining V_k via
%      deterministic least squares.
%
% 3. Matrix Reconstruction:
%    - Apply randomized SVD on the entire first row [X_11, ..., X_1N] to
%      determine row basis vectors.
%    - Use the Hermitian property of X to reconstruct all other rows
%      deterministically.
%
%
%
%
% Inputs:
%   - A: Tensor of size (N, N, M), representing the M initial staes, of size N x N.
%   - b: Measurement result of X, from the observable E_11, (E_12+E_21, ..., E_1N+E_N1), (iE_12-iE_21, ..., iE_1N-iE_N1).
%   - r: Target rank for reconstruction
%   - varargin: Optional arguments (debugON, plotON, X_true, etc, see below.)
%
% Outputs:
%   - X_est: Reconstructed reshaped operator X in matrix form
%   - outputInfo: Structure containing reconstruction details (errors, time, etc.)
%
%
%
% copyright - Quanjun Lang, 2024




%% Parse inputs using inputParser
p = inputParser;

% Required inputs
addRequired(p, 'A');
addRequired(p, 'b');
addRequired(p, 'r');

% Optional arguments with default values
addOptional(p, 'num_retry', 0);         % Number of retries when ALS does not converge
addOptional(p, 'sub_ind_ratio', 0.1)   % ratio of the sub index

addOptional(p, 'plotON', 0);
addOptional(p, 'debugON', 0);
addOptional(p, 'ALS_debugON', 0);
addOptional(p, 'displayON', 0);
addOptional(p, 'compute_dnorm', 0);
addOptional(p, 'Nesterov_beta', 0.5);
addOptional(p, 'ALS_rel_err_tol', 1e-8);
addOptional(p, 'ALS_loss_tol', 1e-8);


addOptional(p, 'X_true', 0);            % True information to compare error
addOptional(p, 'X_true_sub_blocks', 0);


% Parse arguments
parse(p, A, b, r, varargin{:});
% Extract parsed inputs
plotON              = p.Results.plotON;
debugON             = p.Results.debugON;
ALS_debugON         = p.Results.ALS_debugON;
displayON           = p.Results.displayON;
compute_dnorm       = p.Results.compute_dnorm;

X_true_sub_blocks   = p.Results.X_true_sub_blocks;
X_true              = p.Results.X_true;
num_retry           = p.Results.num_retry;
sub_ind_ratio       = p.Results.sub_ind_ratio;

Nesterov_beta        = p.Results.Nesterov_beta;
ALS_rel_err_tol     = p.Results.ALS_rel_err_tol;
ALS_loss_tol        = p.Results.ALS_loss_tol;


%% Extract parameters from input matrices



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


%% Reconstruction of a subset of the first row
% sub_ind_ratio = 0.1;
sub_ind = [1; 1+randsample(N-1, floor((N-1)*sub_ind_ratio))];

% fprintf('Joint ALS for the following blocks in the first row : ');
fprintf('\nLearning operator using ALS on a length = %d random subset of the first row: ', length(sub_ind))
fprintf('%d, ', sort(sub_ind))
fprintf('\n')


small_new_b = new_b(sub_ind, :);
try
    small_X_true = [X_true_sub_blocks{1, sub_ind}];
catch ME
    small_X_true = 0;
end


% ALS_debugON = 1;

[X_est_first_row_small, ALS_outputInfo] = ALS_shared_U(A, small_new_b, r, 'NesterovON', 1, 'X_true', small_X_true, 'debugON', ALS_debugON, 'Nesterov_beta', Nesterov_beta, ...
    'loss_tol', ALS_loss_tol, 'rel_err_tol', ALS_rel_err_tol);

flag = ALS_outputInfo.convergence_flag;
if flag ~= 1 && num_retry>0
    for i = 1:num_retry
        fprintf('%d-th retry: ', i);
        [X_est_first_row_small, ALS_outputInfo] = ALS_shared_U(A, small_new_b, r, 'NesterovON', 1, 'X_true', small_X_true, 'debugON', ALS_debugON, 'Nesterov_beta', Nesterov_beta, ...
            'loss_tol', ALS_loss_tol, 'rel_err_tol', ALS_rel_err_tol);
        flag = ALS_outputInfo.convergence_flag;
        if flag == 1
            break
        end
    end
end


if flag ~= 1
    fprintf('Recovery failed. Please consider generate more data, or increase the retry numbers\n');
    X_est = [];
    outputInfo.rel_error_X = 100;
    outputInfo.flag = flag;
        outputInfo.time = ALS_outputInfo.time;
        outputInfo.ALS_info = ALS_outputInfo;
    return
end

for k = 1:length(sub_ind)
    X_est_blocks{1, sub_ind(k)} = X_est_first_row_small(:, (k-1)*N+1:k*N);
end

%% Reconstruction of the first row
U1 = ALS_outputInfo.all_U(:, :, end);

regmat_V = zeros(M, N*r);
for m = 1:M
    regmat_V(m, :) = vec(A(:, :, m)' * U1)'; % Matrix product A_i * V
end


for k = 1:N
    if ~ any(sub_ind == k)
        b_k = new_b(k, :).';
        vec_V1 = lsqminnorm(regmat_V, conj(b_k));
        Vk = reshape(vec_V1, [N, r]);
        X_est_blocks{1, k} = U1*Vk';
    end
end

X_est_first_row = [X_est_blocks{1, :}];
%% Deterministic reconstruction of the large matrix
X_est_rows = cell(N, 1);
X_est_rows{1} = X_est_first_row;

% [~, ~, V] = svd(X_est_first_row);
% V0 = V(:, 1:r);
% V0_11 = V(1:N, 1:r);
% 
% for k = 2:N
%     X_k1 = X_est_blocks{1, k}';
%     C = (pinv(V0_11)*(X_k1'))';
%     X_est_rows{k} = C * V0';
% end

[~, ~, V0] = rsvd(X_est_first_row, r);   % randomized SVD

V0_11 = V0(1:N, :);
V0_11_inv = pinv(V0_11);

for k = 2:N
    X_k1 = X_est_blocks{1, k}';
    C = (V0_11_inv*(X_k1'))';
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




%% Store estimation results
outputInfo.time = toc(time_counter);
outputInfo.ALS_info = ALS_outputInfo;
outputInfo.flag = flag;


outputInfo.X_est_blocks = X_est_blocks;
outputInfo.X_est_first_row_blocks = X_est_first_row;

fprintf('Finished in %.2f seconds \n', outputInfo.time)
%% Error analysis


outputInfo = error_analysis(outputInfo, X_est, X_true, 'X_est_blocks', X_est_blocks, 'X_true_sub_blocks', X_true_sub_blocks, 'displayON', displayON, 'compute_dnorm', compute_dnorm);





% if ismatrix(X_true)
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
% 
%     outputInfo.error_fro_blocks = error_fro_blocks;
%     outputInfo.rel_error_fro_blocks = rel_error_fro_blocks;
%     outputInfo.norm_blocks = norm_blocks;
% 
% 
% 
%     % Frobenius errors for the entire matrix
%     error_fro_X = norm(X_est - X_true, 'fro');
%     rel_error_fro_X = error_fro_X / norm(X_true, 'fro');
% 
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
%     % % Diamond norm error of Choi matrix
%     % Operator_est = rearrangement_R(X_est);
%     % Choi_est = RR(Operator_est);
%     % 
%     % Operator_true = rearrangement_R(X_true);
%     % Choi_true = RR(Operator_true);
%     % 
%     % outputInfo.error_diamond = dnorm(Choi_est - Choi_true);
%     % 
%     % 
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

% if debugON && ismatrix(X_true)
%     fprintf('Relative Frobenius error \t= %.8f, \tRelative Nuclear error on Choi \t= %.8f, \tRelative Nuclear error on operator \t= %.8f \n',rel_error_fro_X, rel_error_nuc_X, rel_error_nuc_op);
%     fprintf('Frobenius error \t\t= %.8f, \tNuclear error on Choi \t\t= %.8f, \tNuclear error on operator \t\t= %.8f \n',error_fro_X, error_nuc_X, error_nuc_op);
% end

end
