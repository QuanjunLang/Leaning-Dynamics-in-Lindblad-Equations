clc
close all
clear all
rng(0)

addPaths


% copyright - Quanjun Lang, 2024
%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 16;           %
sysInfo.M               = 200;          % number of independent trajectories
sysInfo.dt              = 0.000001;     % true data generation time grid
sysInfo.p               = 1;            % number of jumpoperators
sysInfo.steps           = 3;
sysInfo.channel_dt_rate = 3;




sysInfo.observable_option  = 'First_row_col';
% sysInfo.observable_option  = 'First_row_col_diag';
% sysInfo.observable_option  = 'Full_state';
% sysInfo.observable_option  = 'Multiple_random_observables';  sysInfo.N_o = 7;
% sysInfo.observable_option  = 'Single_random_observable';

sysInfo = update_sys(sysInfo);

sysInfo.PAPER_FIG_DIR = 'figure';
%%
[all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);


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



%% example: one block of L (high rank Hermittian)
ind_o = 1;
b_L = squeeze(all_rho(ind_o, 2, :) - all_rho(ind_o, 1, :))/sysInfo.dt;
% b_L_est = squeeze(sum(conj(RL_obs_blocks{ind_o}).*all_rho0, [1,2]));
l = RL_obs_blocks{ind_o};

[L_blocks_ind_o_est, outputInfo_L] = ALS(all_rho0, b_L, rank(RL_obs_blocks{ind_o}), 'X_true', RL_obs_blocks{ind_o}, 'debugON', 1, 'NesterovON', 1);
[L_blocks_ind_o_est, outputInfo_L] = ALS(all_rho0, b_L, rank(RL_obs_blocks{ind_o}), 'X_true', RL_obs_blocks{ind_o}, 'debugON', 1, 'NesterovON', 0);


%% example: one block of L (low rank non-Hermittian)
ind_first_row = 1;
l = RL_sub_blocks{1, ind_first_row+1};
b = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
b_fisrt_row = (b(n+ind_first_row, :) -1i*b(2*n-1+ind_first_row, :)).'/2;
[L_blocks_ind_o_est, outputInfo_L] = ALS(all_rho0, b_fisrt_row, rank(l), 'X_true', l, 'debugON', 1, 'NesterovON', 0);
[L_blocks_ind_o_est, outputInfo_L] = ALS(all_rho0, b_fisrt_row, rank(l), 'X_true', l, 'debugON', 1, 'NesterovON', 1);
%% example: one block of E
% ind_o = 10;
%
% b_E = squeeze(all_rho(ind_o, sysInfo.channel_dt_rate+1, :));
%
%
% [E_blocks_ind_o_est, ~] = ALS(all_rho0, b_E, rank(RE_obs_blocks{ind_o}), 'X_true', RL_obs_blocks{ind_o}, 'debugON', 1);

%% Parallel ALS for L, first only
A = all_rho0;
bL = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
rL = trueInfo.rank_RL_true;

% [RL_est, RL_Info] = parallel_ALS_first_row(A, bL, rL, 'X_true_sub_blocks', RL_sub_blocks, 'X_true_obs_blocks', RL_obs_blocks, 'X_true', RL, 'num_retry', 0);
[RL_est, RL_Info] = parallel_ALS_first_row_new(A, bL, rL, 'X_true_sub_blocks', RL_sub_blocks, 'X_true_obs_blocks', RL_obs_blocks, 'X_true', RL, 'num_retry', 0);

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


    fprintf('Parallel learning L, error = %.8f, time = %.2f seconds \n', RL_Info.rel_error_X, RL_Info.time)
end

%% Parallel ALS for L, first row and diagonal
% A = all_rho0;
% bL = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
% rL = trueInfo.rank_RL_true;
%
% [RL_est, RL_Info] = parallel_ALS_low_rank(A, bL, rL, 'X_true_sub_blocks', RL_sub_blocks, 'X_true_obs_blocks', RL_obs_blocks, 'X_true', RL);
%
% figure;
% subplot(131)
% surf(RL_Info.error_blocks);view(0, 90);
% colorbar
% title('L blocks error')
% subplot(132)
% surf(RL_Info.rel_error_blocks);view(0, 90);
% colorbar
% title('L blocks relative error')
% subplot(133)
% surf(RL_Info.norm_blocks);view(0, 90);
% colorbar
% title('L blocks norm')
%
%
% fprintf('Parallel learning L, error = %.8f, time = %.2f seconds \n', RL_Info.rel_error_X, RL_Info.time)

%% Parallel ALS first row, first column, diagonal
% A = all_rho0;
% bL = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
% rL = trueInfo.rank_RL_true;
%
% [RL_est, RL_Info] = parallel_ALS(A, bL, rL, 'X_true_sub_blocks', RL_sub_blocks, 'X_true_obs_blocks', RL_obs_blocks, 'X_true', RL);
%
% figure;
% subplot(131)
% surf(RL_Info.error_blocks);view(0, 90);
% colorbar
% title('L blocks error')
% subplot(132)
% surf(RL_Info.rel_error_blocks);view(0, 90);
% colorbar
% title('L blocks relative error')
% subplot(133)
% surf(RL_Info.norm_blocks);view(0, 90);
% colorbar
% title('L blocks norm')
%
%
% fprintf('Parallel learning L, error = %.8f, time = %.2f seconds \n', RL_Info.rel_error_X, RL_Info.time)

%% Parallel ALS for E
% A = all_rho0;
% bE = squeeze(all_rho(:, sysInfo.channel_dt_rate+1, :));
% rE = trueInfo.rank_RE_true;
%
% [RE_est, RE_Info] = parallel_ALS(A, bE, rE, 'X_true_sub_blocks', RE_sub_blocks);
%
% norm(RE_est - RE, 'fro')
%
% figure;
% surf(RE_Info.error_blocks);view(0, 90);
% colorbar
% title('E blocks error')

%%

% all_L_blocks_svd = zeros(n, N_o);
% all_L_blocks_rank = zeros(N_o, 1);
% for ind_o = 1:N_o
%     all_L_blocks_svd(:, ind_o) = svd(RL_obs_blocks{ind_o});
%     all_L_blocks_rank(ind_o) = rank(RL_obs_blocks{ind_o});
% end
%
% all_E_blocks_svd = zeros(n, N_o);
% all_E_blocks_rank = zeros(N_o, 1);
% for ind_o = 1:N_o
%     all_E_blocks_svd(:, ind_o) = svd(RE_obs_blocks{ind_o});
%     all_E_blocks_rank(ind_o) = rank(RE_obs_blocks{ind_o});
% end
%
%
% figure;
% subplot(121);hold on;
% plot(all_L_blocks_rank, 'b')
% plot(all_E_blocks_rank, 'r')
% title('rank of subblocks of E and L')
% legend('L', 'E')
% xlabel('index of observables')
% subplot(122);hold on;
% plot(log10(all_L_blocks_svd), 'b')
% plot(log10(all_E_blocks_svd), 'r')
% title('log10 svd of subblocks of E and L')
% xlabel('index of observables')


