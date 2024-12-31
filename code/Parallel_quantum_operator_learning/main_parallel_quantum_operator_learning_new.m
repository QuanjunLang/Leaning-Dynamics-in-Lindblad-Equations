clc
close all
clear all
rng(0)

addPaths


% copyright - Quanjun Lang, 2024
%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 32;           %
sysInfo.M               = 800;          % number of independent trajectories
sysInfo.dt              = 0.0000001;     % true data generation time grid
sysInfo.p               = 2;            % number of jumpoperators
sysInfo.steps           = 3;
sysInfo.channel_dt_rate = 3;



% sysInfo.observable_option  = 'Channel_first_row_col';
sysInfo.observable_option  = 'First_row_col';
% sysInfo.observable_option  = 'First_row_col_diag';
% sysInfo.observable_option  = 'Full_state';
% sysInfo.observable_option  = 'Multiple_random_observables';  sysInfo.N_o = 7;
% sysInfo.observable_option  = 'Single_random_observable';

sysInfo = update_sys(sysInfo);

sysInfo.PAPER_FIG_DIR = 'figure';
%%
if contains(sysInfo.observable_option, 'Channel')
    [all_rho, trueInfo, observableInfo] = generate_data_channel(sysInfo, 'plotON', 1);
else
    [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);
end


%% plot rank of the chennel operator with different time steps

all_dt = 10.^linspace(-7, -1, 7);
s = zeros(length(all_dt), sysInfo.n^2);
r = zeros(length(all_dt), 1);
for i = 1:length(all_dt)
    dt = all_dt(i);
    E_temp = rearrangement_R(expm(trueInfo.L_true * dt));
    r(i) = rank(E_temp);
    s(i, :) = svd(E_temp);
end

figure;hold on;
set(gca,'DefaultLineLineWidth',2)
plot(log10(s'))
xline(trueInfo.rank_RL_true, 'r-')
xline(trueInfo.rank_RE_true, 'b-')
xline(sysInfo.n, 'k--')

for i = 1:length(all_dt)
    plot(r(i), log10(s(i, r(i))), 'o');
end

legendStrings = arrayfun(@num2str, all_dt, 'UniformOutput', false);
legend([legendStrings, 'rank(RL)', ['rank(RE), dt = ', num2str(sysInfo.channel_dt)], 'n'], 'Location', 'best');

xlabel('index of singular values')
ylabel('log_{10} of singular values')
title('The change of channel operator rank with different time scale')

%% Multiple block ALS
n = sysInfo.n;
M = sysInfo.M;
N_o = sysInfo.N_o;


L = trueInfo.L_true;
E = trueInfo.E_true;
RL = trueInfo.RL_true;
RE = trueInfo.RE_true;
rank_RL = trueInfo.rank_RL_true;
rank_RE = trueInfo.rank_RE_true;

RE_sub_blocks = trueInfo.RE_sub_blocks;
RL_sub_blocks = trueInfo.RL_sub_blocks;
RE_obs_blocks = trueInfo.RE_obs_blocks;
RL_obs_blocks = trueInfo.RL_obs_blocks;

all_rho0    = observableInfo.rho0;


%% example: one block of L (low rank non-Hermittian)
b = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
new_b = zeros(n, M);
new_b(1, :) = b(1, :);
for i = 1:n-1
    new_b(i+1, :) = (b(1+i, :) -1i*b(n+i, :))/2;
end

ind_first_row = 1;
l = RL_sub_blocks{ind_first_row, 1};

ALS(all_rho0, new_b(ind_first_row, :).', rank(l), 'X_true', l, 'debugON', 1, 'NesterovON', 0);
ALS(all_rho0, new_b(ind_first_row, :).', rank(l), 'X_true', l, 'debugON', 1, 'NesterovON', 1);

%% example: learning the first row of L
[X_est_first_row, outputInfo] = ALS_shared_U(all_rho0, new_b, rank_RL, 'NesterovON', 1, 'X_true', trueInfo.RL_true(1:n, :), 'debugON', 1);



%% Parallel ALS for L, first only
A = all_rho0;
bL = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;


[RL_est, RL_Info] = parallel_ALS_first_row_new(A, bL, rank_RL, 'X_true_sub_blocks', RL_sub_blocks, 'X_true_obs_blocks', RL_obs_blocks, 'X_true', RL, 'num_retry', 0);

if RL_Info.flag
    figure;
    subplot(131)
    % surf(RL_Info.error_blocks);view(0, 90);
    imagesc(RL_Info.error_blocks)
    colorbar
    title('L blocks error')
    subplot(132)
    % surf(RL_Info.rel_error_blocks);view(0, 90);
    imagesc(RL_Info.rel_error_blocks)
    colorbar
    title('L blocks relative error')
    subplot(133)
    % surf(RL_Info.norm_blocks);view(0, 90);
    imagesc(RL_Info.norm_blocks)
    colorbar
    title('L blocks norm')
    

    % fprintf('Parallel learning L, relative Frobenius error = %.8f, relative Nuclear error = %.8f, time = %.2f seconds \n', RL_Info.rel_error_X, RL_Info.rel_error_X_nuclear, RL_Info.time)
    fprintf('Parallel learning L, relative Frobenius error = %.8f, relative Nuclear error on Choi = %.8f, Frobenius error = %.8f, Nuclear error on Choi = %.8f, time = %.2f seconds \n' ...
        , RL_Info.rel_error_X, RL_Info.rel_error_X_nuclear, RL_Info.error_X, RL_Info.error_X_nuclear, RL_Info.time)
    fprintf('Nuclear norm on operator = %.8f, relative Nuclear norm on operator=%.8f', RL_Info.error_X_nuclear_op, RL_Info.rel_error_X_nuclear_op)

end




%% Learning L only from blocks on the first row

A = all_rho0;
bL = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;


[RL_est_first_row, RL_est_first_row_Info] = ALS_first_row_joint_learning(A, b, rank_RL, 'X_true_sub_blocks', RL_sub_blocks, 'X_true', RL, 'num_retry', 0);



if RL_est_first_row_Info.flag
    figure;
    subplot(131)
    % surf(RL_Info.error_blocks);view(0, 90);
    imagesc(RL_est_first_row_Info.error_blocks)
    colorbar
    title('L blocks error')
    subplot(132)
    % surf(RL_Info.rel_error_blocks);view(0, 90);
    imagesc(RL_est_first_row_Info.rel_error_blocks)
    colorbar
    title('L blocks relative error')
    subplot(133)
    % surf(RL_Info.norm_blocks);view(0, 90);
    imagesc(RL_est_first_row_Info.norm_blocks)
    colorbar
    title('L blocks norm')


    % fprintf('Parallel learning L, relative Frobenius error = %.8f, relative Nuclear error = %.8f, time = %.2f seconds \n', RL_est_first_row_Info.rel_error_X, RL_est_first_row_Info.rel_error_X_nuclear, RL_est_first_row_Info.time)


    fprintf('First row joint learning, relative Frobenius error = %.8f, relative Nuclear error on Choi = %.8f, Frobenius error = %.8f, Nuclear error on Choi = %.8f, time = %.2f seconds \n' ...
        , RL_est_first_row_Info.rel_error_X, RL_est_first_row_Info.rel_error_X_nuclear, RL_est_first_row_Info.error_X, RL_est_first_row_Info.error_X_nuclear, RL_est_first_row_Info.time)
    fprintf('Nuclear norm on operator = %.8f, relative Nuclear norm on operator=%.8f', RL_est_first_row_Info.error_X_nuclear_op, RL_est_first_row_Info.rel_error_X_nuclear_op)


end


%%

% Nuclear norm error of operator
error_X_nuclear_op     = norm(svd(rearrangement_R(RL_est_first_row) - rearrangement_R(RL)), 1)
rel_error_X_nuclear_op = error_X_nuclear_op./sum(svd(rearrangement_R(RL)))
% 
% outputInfo.error_X_nuclear_op = error_X_nuclear_op;
% outputInfo.rel_error_X_nuclear_op = rel_error_X_nuclear_op;

%%
% RL_est_first_row_blocks     = RL_Info.X_est_first_row_blocks;
% RL_est_sub_blocks           = RL_Info.X_est_blocks;
% RL_true_first_row_blocks    = trueInfo.RL_true(1:sysInfo.n, :);
% 
% sum(svd(RL_est_first_row_blocks - RL_true_first_row_blocks))
% 
% sum(svd(RL_est - trueInfo.RL_true))
% 
% 
% %%
% b = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
% new_b = zeros(n, M);
% new_b(1, :) = b(1, :);
% for i = 1:n-1
%     new_b(i+1, :) = (b(1+i, :) -1i*b(n+i, :))/2;
% end
% 
% [U, S, V] = svd(RL_est_first_row_blocks);
% 
% U0 = U(:, 1:r);
% V0 = (S(1:r, 1:r)*V(:, 1:r)')';
% 
% tic
% % [M_est, outputInfo] = ALS_shared_U(A, new_b, rL, 'U0', U0, 'V0', V0, 'NesterovON', 0, 'maxIter', 10);
% [M_est, outputInfo] = ALS_shared_U(A, new_b, rL, 'NesterovON', 1);
% toc
% 
% 
% 
% %%
% r = trueInfo.rank_RL_true;
% [U, S, V] = svd(M_est);
% U0 = U(:, 1:r);
% S0 = S(1:r, 1:r);
% V0 = V(:, 1:r);
% V0_11 = V(1:n, 1:r);
% X_11 = M_est(1:n, 1:n);% Current off-diagonal block
% % Reconstruct the remaining blocks using SVD and basis projections
% 
% 
% 
% 
% for k = 2:n
%     X_k1 = M_est(:, (k-1)*n+1:k*n)'; % Diagonal block for the k-th row/column
%     C = (pinv(V0_11)*(X_k1'))';
% 
%     M_est_blocks{k} = C * V0';
% 
% end
% 
% 




%%
% figure;
% plot(log10(outputInfo.loss))
% 
% 
% sum(svd(M_est - RL_true_first_row_blocks))
% 
