clc
close all
clear all
rng(0)
addPaths




% copyright - Quanjun Lang, 2024
%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 8;           %
sysInfo.M               = 50;          % number of independent trajectories
sysInfo.p               = 3;           % number of jumpoperators


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



sigma = 1e-4;



Num_samples = 10;



all_ALS_N2 = cell(Num_samples, 1);
all_ALS_P = cell(Num_samples, 1);
all_ALS_N = cell(Num_samples, 1);
all_ALS_I = cell(Num_samples, 1);


for i = 1:Num_samples
    if contains(sysInfo.observable_option, 'Channel')
        [all_rho, trueInfo, observableInfo] = generate_data_channel(sysInfo, 'plotON', 1);
    else
        [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);
    end

    
    all_rho = all_rho + randn(size(all_rho))*sigma;

    %% Multiple block ALS
    n = sysInfo.n;
    M = sysInfo.M;
    N_o = sysInfo.N_o;

    E = trueInfo.E_true;
    RE = trueInfo.RE_true;
    rank_RE = trueInfo.rank_RE_true;
    RE_sub_blocks = trueInfo.RE_sub_blocks;
    all_rho0    = observableInfo.rho0;



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

    fprintf('\nJoint ALS for the entire matrix...\n')
    [RE_est_total, RE_Info_total] = ALS(A_mat, b_mat, rank_RE, 'X_true', RE, 'debugON', 0, 'Nesterov_beta', 0.2, 'rel_err_tol', 1e-5);
    fprintf('Finished in %.2f seconds \n', RE_Info_total.time)

    RE_Info_total = error_analysis(RE_Info_total, RE_est_total, RE, 'displayON', 1, 'compute_dnorm', 1);

    %% Learning in parallel for all first row blocks
    [RE_est, RE_Info] = parallel_ALS_first_row_new(all_rho0, all_rho, rank_RE, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 0, 'displayON', 1, 'compute_dnorm', 1, 'Nesterov_beta', 0.4, 'num_retry', 3);


    %% Learning the first row jointly
    [RE_est_first_row, RE_est_first_row_Info] = ALS_first_row_joint_learning(all_rho0, all_rho, rank_RE, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'num_retry', 3, 'debugON', 0, 'displayON', 1, 'compute_dnorm', 1, 'Nesterov_beta', 0.5);


    %% Learning a random subset of the first row jointly
    [RE_est_first_row_subset, RE_est_first_row_Info_subset] = ALS_first_row_joint_learning_random_subset(all_rho0, all_rho, rank_RE, 'sub_ind_ratio', 0.2, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 0, 'plotON', 0, 'ALS_debugON', 0, 'displayON', 1, 'compute_dnorm', 1);


    %% store the numbers
    all_ALS_N2{i} = RE_Info_total;
    all_ALS_P{i} = RE_Info;
    all_ALS_N{i} = RE_est_first_row_Info;
    all_ALS_I{i} = RE_est_first_row_Info_subset;

end


%% prepare the numbers

Num_items = 3;

all_error = zeros(4, Num_items, Num_samples);

for i = 1:Num_samples
    all_error(1, 1, i) = all_ALS_N2{i}.rel_error_fro_X;
    all_error(1, 2, i) = all_ALS_N2{i}.error_diamond;
    all_error(1, 3, i) = all_ALS_N2{i}.time;

    all_error(2, 1, i) = all_ALS_P{i}.rel_error_fro_X;
    all_error(2, 2, i) = all_ALS_P{i}.error_diamond;
    all_error(2, 3, i) = all_ALS_P{i}.time;

    all_error(3, 1, i) = all_ALS_N{i}.rel_error_fro_X;
    all_error(3, 2, i) = all_ALS_N{i}.error_diamond;
    all_error(3, 3, i) = all_ALS_N{i}.time;

    all_error(4, 1, i) = all_ALS_I{i}.rel_error_fro_X;
    all_error(4, 2, i) = all_ALS_I{i}.error_diamond;
    all_error(4, 3, i) = all_ALS_I{i}.time;
end

%%

methods = {'ALS-$N^2$', 'ALS-P', 'ALS-N', 'ALS-I'};
mean1 = mean(squeeze(all_error(:, 1, :))');
std1 = std(squeeze(all_error(:, 1, :))');
mean2 = mean(squeeze(all_error(:, 2, :))');
std2 = std(squeeze(all_error(:, 2, :))');
mean3 = mean(squeeze(all_error(:, 3, :))');
std3 = std(squeeze(all_error(:, 3, :))');


frob_err = strcat(num2str(mean1', '%.2e'), ' $\pm$ ', num2str(std1', '%.2e'));
diamond_err = strcat(num2str(mean2', '%.2e'), ' $\pm$ ', num2str(std2', '%.2e'));
time_err = strcat(num2str(mean3', '%.2e'), ' $\pm$ ', num2str(std3', '%.2e'));


% Combine into a table
T = table(methods', frob_err, diamond_err, time_err, ...
    'VariableNames', {'Method', 'Frobenius Error', 'Diamond Error', 'Time(s)'});

% Write to CSV
writetable(T, '5_1_1_channel_typical_data.csv')
writetable(T, '/Users/quanjunlang/Documents/GitHub/Learning-operators-in-Lindblad-Master-Equation-Paper/5_1_1_channel_typical_data.csv')


% currentFile = mfilename('fullpath');
% 
% data = table(methods', mean1', std1', mean2', std2', mean3', std3', ...
%     'VariableNames', {'Method', 'Mean1', 'Std1', 'Mean2', 'Std2', 'Mean3', 'Std3'});
% writetable(data, 'data.csv');


% %%
% if contains(sysInfo.observable_option, 'Channel')
%     [all_rho, trueInfo, observableInfo] = generate_data_channel(sysInfo, 'plotON', 1);
% else
%     [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);
% end
%
% sigma = 1e-4;
% all_rho = all_rho + randn(size(all_rho))*sigma;
%
% %% Multiple block ALS
% n = sysInfo.n;
% M = sysInfo.M;
% N_o = sysInfo.N_o;
%
% E = trueInfo.E_true;
% RE = trueInfo.RE_true;
% rank_RE = trueInfo.rank_RE_true;
% RE_sub_blocks = trueInfo.RE_sub_blocks;
% all_rho0    = observableInfo.rho0;
%
%
%

% %% Total ALS
%
% % single observable case
% % A_mat = zeros(N^2, N^2, M);
% % for m = 1:M
% %     O = all_O(:, :, m);
% %     rho0 = all_rho0(:, :, m);
% %     A_mat(:, :, m) = kron(conj(O), rho0);
% % end
%
% % multiple shared observables
% A_mat = zeros(n^2, n^2, M*N_o);
% l = 1;
% for m = 1:M
%     rho0 = all_rho0(:, :, m);
%     for k = 1:N_o
%         O = observableInfo.O(:, :, k);
%         A_mat(:, :, l) = kron(conj(O), rho0);
%         l = l+1;
%     end
% end
%
% b_mat = vec(all_rho);
%
% fprintf('\nJoint ALS for the entire matrix...\n')
% [RE_est_total, RE_Info_total] = ALS(A_mat, b_mat, rank_RE, 'X_true', RE, 'debugON', 0);
% fprintf('Finished in %.2f seconds \n', RE_Info_total.time)
%
% RE_Info_total = error_analysis(RE_Info_total, RE_est_total, RE, 'displayON', 1, 'compute_dnorm', 1);
%
% %% Learning in parallel for all first row blocks
% [RE_est, RE_Info] = parallel_ALS_first_row_new(all_rho0, all_rho, rank_RE, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 0, 'displayON', 1, 'compute_dnorm', 1);
%
%
% %% Learning the first row jointly
% [RE_est_first_row, RE_est_first_row_Info] = ALS_first_row_joint_learning(all_rho0, all_rho, rank_RE, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'num_retry', 0, 'debugON', 0, 'displayON', 1, 'compute_dnorm', 1);
%
%
% %% Learning a random subset of the first row jointly
% [RE_est_first_row_subset, RE_est_first_row_Info_subset] = ALS_first_row_joint_learning_random_subset(all_rho0, all_rho, rank_RE, 'sub_ind_ratio', 0.2, 'X_true_sub_blocks', RE_sub_blocks, 'X_true', RE, 'debugON', 0, 'plotON', 0, 'ALS_debugON', 0, 'displayON', 1, 'compute_dnorm', 1);
%
