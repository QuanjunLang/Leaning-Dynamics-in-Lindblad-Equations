clc
close all
clear all
rng(0)
addPaths

% copyright - Quanjun Lang, 2024
%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 12;           %
sysInfo.M               = 50;          % number of independent trajectories
sysInfo.p               = 2;            % number of jumpoperators


sysInfo.steps = 1;
sysInfo.dt = 1;
sysInfo.channel_dt_rate = 1;

sysInfo.observable_option  = 'Channel_first_row_col';
% sysInfo.observable_option  = 'First_row_col';
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

sigma = 1e-5;
all_rho = all_rho + randn(size(all_rho))*sigma;

%% Multiple block ALS
n = sysInfo.n;
M = sysInfo.M;
N_o = sysInfo.N_o;

E = trueInfo.E_true;
RE = trueInfo.RE_true;
rank_RE = trueInfo.rank_RE_true;
RE_sub_blocks = trueInfo.RE_sub_blocks;
% RE_obs_blocks = trueInfo.RE_obs_blocks;
all_rho0    = observableInfo.rho0;

%% example: one block of L (low rank non-Hermittian)
% b = all_rho;
% new_b = zeros(n, M);
% new_b(1, :) = b(1, :);
% for i = 1:n-1
%     new_b(i+1, :) = (b(1+i, :) -1i*b(n+i, :))/2;
% end
% % 
% % ind_first_row = 1;
% % l = RE_sub_blocks{ind_first_row, 1};
% % 
% % ALS(all_rho0, new_b(ind_first_row, :).', rank(l), 'X_true', l, 'debugON', 1, 'NesterovON', 0);
% % ALS(all_rho0, new_b(ind_first_row, :).', rank(l), 'X_true', l, 'debugON', 1, 'NesterovON', 1);

%% example: learning the first row

% sub_ind = randsample(n, floor(log(n)));
% 
% small_new_b = new_b(sub_ind, :);
% small_X_true = [RE_sub_blocks{1, sub_ind}];
% [X_est_first_row_small, outputInfo_first_row_small] = ALS_shared_U(all_rho0, small_new_b, rank_RE, 'NesterovON', 1, 'X_true', small_X_true, 'debugON', 1);
% 
% 
% 
% [X_est_first_row, outputInfo] = ALS_shared_U(all_rho0, new_b, rank_RE, 'NesterovON', 1, 'X_true', trueInfo.RE_true(1:n, :), 'debugON', 1);


% [X_1, outputInfo_1] = ALS(all_rho0, small_new_b.', rank_RE, 'NesterovON', 1, 'X_true', small_X_true, 'debugON', 1);


%% Total ALS 

% single observable case
% A_mat = zeros(N^2, N^2, M);
% for m = 1:M
%     O = all_O(:, :, m);
%     rho0 = all_rho0(:, :, m);
%     A_mat(:, :, m) = kron(conj(O), rho0);
% end

% multiple shared observables
A_mat = zeros(n^2, n^2, M*N_o);
l = 1;
for m = 1:M
    rho0 = all_rho0(:, :, m);
    for k = 1:N_o
        O = observableInfo.O(:, :, k);
        A_mat(:, :, l) = kron(conj(O), rho0);
        l = l+1;
    end
end

b_mat = vec(all_rho);

[RE_est_total, RE_Info_total] = ALS(A_mat, b_mat, rank_RE, 'X_true', RE, 'debugON', 1);
fprintf('Finished in %.2f seconds \n', RE_Info_total.time)

RE_Info_total = error_analysis(RE_Info_total, RE_est_total, RE, 'displayON', 1, 'compute_dnorm', 0);

%% Learning in parallel for all first row blocks
[RE_est, RE_Info] = parallel_ALS_first_row_new(all_rho0, all_rho, rank_RE, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 1, 'displayON', 1, 'compute_dnorm', 0);


%% Learning the first row jointly
[RE_est_first_row, RE_est_first_row_Info] = ALS_first_row_joint_learning(all_rho0, all_rho, rank_RE, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'num_retry', 0, 'debugON', 1, 'displayON', 1, 'compute_dnorm', 0);


%% Learning a random subset of the first row jointly
[RE_est_first_row_subset, RE_est_first_row_Info_subset] = ALS_first_row_joint_learning_random_subset(all_rho0, all_rho, rank_RE, 'sub_ind_ratio', 0.4, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 1, 'plotON', 0, 'ALS_debugON', 0, 'displayON', 1, 'compute_dnorm', 0);



%%
K = sysInfo.n/5;
L = 3;

all_ratio = linspace(0, 1, K);
all_error_fro = zeros(K, L);
all_error_nuc = zeros(K, L);


all_time = zeros(K, L);
all_iters = zeros(K, L);
all_time_ALS = zeros(K, L);


all_flags = zeros(K, L);


for i = 1:K
    for l = 1:L
        sub_ind_ratio = all_ratio(i);
        [RE_est_first_row_subset, RE_est_first_row_Info_subset] = ALS_first_row_joint_learning_random_subset(all_rho0, all_rho, rank_RE, 'sub_ind_ratio', sub_ind_ratio, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 0, 'plotON', 0, 'ALS_debugON', 0);
        if RE_est_first_row_Info_subset.flag
            all_error_fro(i, l) = RE_est_first_row_Info_subset.rel_error_fro_X;
            all_error_nuc(i, l) = RE_est_first_row_Info_subset.rel_error_nuc_X;
            all_time(i, l) = RE_est_first_row_Info_subset.time;
            all_iters(i, l) = RE_est_first_row_Info_subset.ALS_info.iters;
            all_time_ALS(i, l) = RE_est_first_row_Info_subset.ALS_info.time;
            all_flags(i, l) = 1;
        else
            all_error_fro(i, l) = 10^3;
            all_error_nuc(i, l) = 10^3;
            all_time(i, l) = 10^3;
            all_time_ALS(i, l) = 10^3;
            all_iters(i, l) = 10^3;
            all_flags(i, l) = 0;
        end
    end
end



%%
FLAG = prod(all_flags, 2);
ind_temp = 1:K;

conv_ind = ind_temp(FLAG >0);


%%
figure;
subplot(131);hold on;
plot(all_ratio(conv_ind), log10(all_error_fro(conv_ind, :)), 'Color',color1);
plot(all_ratio(conv_ind), log10(all_error_nuc(conv_ind, :)), 'Color',color2);

ymin = min([all_error_fro(conv_ind, :), all_error_nuc(conv_ind, :)],[], 'all');
ymax = max([all_error_fro(conv_ind, :), all_error_nuc(conv_ind, :)],[], 'all');
ylim([log10(ymin*0.9), log10(ymax*1.1)])

legend('Rel Frobenius norm', 'Rel Nuclear norm')

subplot(132);hold on;
plot(all_ratio(conv_ind), all_time(conv_ind, :), 'Color',color1);
plot(all_ratio(conv_ind), all_time_ALS(conv_ind, :), 'Color',color2);

ymin = min([all_time(conv_ind, :), all_time(conv_ind, :)],[], 'all');
ymax = max([all_time_ALS(conv_ind, :), all_time_ALS(conv_ind, :)],[], 'all');
ylim([ymin*0.9, ymax*1.1])

legend('Total time', 'ALS time')


subplot(133)
plot(all_ratio(conv_ind), all_iters(conv_ind, :), 'Color',color1);

ymin = min([all_iters(conv_ind, :)],[], 'all');
ymax = max([all_iters(conv_ind, :)],[], 'all');
ylim([ymin*0.9, ymax*1.1])

legend('Number of iterations')













%%
% dnorm( RE_est_first_row_subset - RE)

%%
% A = randn(4, 4) + 1i*rand(4, 4);
% B = randn(4, 4) + 1i*rand(4, 4);
% C = kron(A, B);
% RC = rearrangement_R(C);
% 
% RC - vec(A.')*vec(B.').'
% 
% 
% RR(C) - vec(B) * vec(A).'
% 
% %% 
% E = trueInfo.E_true;
% RE = rearrangement_R(E);
% Choi = trueInfo.Choi_true;
% 
% RRE = RR(E);
%%
% n = sysInfo.n
% A = randn(4, 4);
% T = zeros(4, 4); T(3, 2) = 1;
% 
% A*T
% 
% %%
% rho0_temp = all_rho0(:, :, 1);
% O_temp = (observableInfo.O(:, :, 4) - 1i*observableInfo.O(:, :, 7))/2;  
% % O_temp = zeros(4, 4);
% % O_temp(2, 1) = 1;
% 
% E_kl = zeros(4, 4);
% E_kl(1, 4) = 1;
% 
% T_temp = kron(O_temp.', rho0_temp');
% sum(conj(RE) .*T_temp, 'all')
% 
% % observableInfo.all_rho1(:, :, 1)
% sum(conj(RE_sub_blocks{4, 1}) .*rho0_temp, 'all')
% sum(RE_sub_blocks{1, 4}.' .*rho0_temp, 'all')
% sum(conj(rho0_temp).*RE_sub_blocks{1, 4}, 'all')
% 
% sum(conj(observableInfo.all_rho1(:, :, 1)).*E_kl, 'all')
% 
% % sum(RE_sub_blocks{1, 4} .*rho0_temp, 'all')
% % observableInfo.all_rho1(:, :, 1) + observableInfo.all_rho1(:, :, 1).'
% % 1i*observableInfo.all_rho1(:, :, 1) -1i* observableInfo.all_rho1(:, :, 1).'
% %%
% rho0_temp = all_rho0(:, :, 1);
% O_temp = observableInfo.O(:, :, 3);
% 
% T_temp = kron(conj(rho0_temp), O_temp);
% sum(conj(T_temp) .*RRE, 'all')
% 
% 
% 
% %%
% x = linspace(40, 10000, 1000);
% x = x.^2;
% c = 10000
% y = (log(x) + c)./(log(log(x)));
% 
% figure;plot(x, y)