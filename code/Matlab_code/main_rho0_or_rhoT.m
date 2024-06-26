clc
close all
clear all
rng(2)
% pyenv(Version='/opt/anaconda3/envs/py3.11/bin/python');
% pyenv(Version='/opt/anaconda3/bin/python');
% pyenv('Version', '');
pyenv(Version='/usr/bin/python3');
% pyenv('ExecutionMode', 'OutOfProcess');

addPaths

%% system settings
sysInfo.n       = 8;            % 
sysInfo.M       = 8;            % number of independent trajectories
sysInfo.dt      = 0.005;        % true data generation time grid
sysInfo.p       = 3;            % number of jump operators

sysInfo.steps   = 1000;
sysInfo = update_sys(sysInfo);

[all_rho, trueInfo] = generate_data(sysInfo);

%% decomposition of the true L into H and K

L_decomp_true = L_decomposition_hamiltonian_kossakowski(trueInfo.L_true, sysInfo.p, trueInfo);
%% generate observational data

obsInfo.obs_std = 1e-5;
obsInfo.obs_gap = 100;
obsInfo.obs_len = 10;


[all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);


%% Demo
plot_sample_traj(sysInfo)


%% Prony fitting test
all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);

ind_1 = 1;
ind_2 = 2;
m = 1;

plot_prony_rho(all_rho_prony, all_rho, all_rho_obs, sysInfo, obsInfo, ind_1, ind_2, m)
% plot_all_prony_modes(all_rho_prony, sysInfo, trueInfo)
%% Learn the channel operators using observed derivatives and prony-fitted derivatives

M = sysInfo.M;
n = sysInfo.n;


% Prony fit data pair
all_rho_prony_pair_data = zeros(n, n, 2, M*(obsInfo.obs_len+1));

for t = 1:obsInfo.obs_len
    all_rho_prony_pair_data(:, :, 1, (t-1)*M+1:t*M) = all_rho_obs(:, :, t, :);
    for i = 1:n
        for j = 1:n
            for m = 1:M
                all_rho_prony_pair_data(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(obsInfo.tgrid(t));
            end
        end
    end
end

% Observation data pair
all_rho_obs_pair_data = zeros(n, n, 2, M*(obsInfo.obs_len+1));
for t = 1:obsInfo.obs_len
    all_rho_obs_pair_data(:, :, 1, (t-1)*M+1:t*M) = all_rho_obs(:, :, t, :);
    all_rho_obs_pair_data(:, :, 2, (t-1)*M+1:t*M) = (all_rho_obs(:, :, t+1, :) - all_rho_obs(:, :, t, :))/obsInfo.dt;
end

% Prony fit data pair augmented
all_rho_prony_pair_data_aug = zeros(n, n, 2, M*(sysInfo.steps+1));

for t = 1:obsInfo.obs_len
    for i = 1:n
        for j = 1:n
            for m = 1:M
                all_rho_prony_pair_data_aug(i, j, 1, (t-1)*M+m) = all_rho_prony{i, j, m}.h(sysInfo.tgrid(t));
                all_rho_prony_pair_data_aug(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(sysInfo.tgrid(t));
            end
        end
    end
end

% Prony fit data time 0
all_rho_prony_pair_data_aug = zeros(n, n, 2, M*(sysInfo.steps+1));

t = 1;
for i = 1:n
    for j = 1:n
        for m = 1:M
            all_rho_prony_pair_t0(i, j, 1, (t-1)*M+m) = all_rho_prony{i, j, m}.h(sysInfo.tgrid(t));
            all_rho_prony_pair_t0(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(sysInfo.tgrid(t));
        end
    end
end


%% Two stage
% get L using Prony fitted data
result_prony = multi_mat_als(all_rho_prony_pair_data, trueInfo.r_true, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);

% get L using finite difference derivative
result_obs = multi_mat_als(all_rho_obs_pair_data, trueInfo.r_true, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);

% get L using time 0 only
result_t0 = multi_mat_als(all_rho_prony_pair_t0, trueInfo.r_true, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);


%%
figure;

hold on; grid on
plot(log10(result_prony.loss), 'DisplayName','Loss Prony');
plot(log10(result_prony.err_E), ':', 'linewidth', 3, 'DisplayName','E Prony');

plot(log10(result_obs.loss), 'DisplayName','Loss obs');
plot(log10(result_obs.err_E), ':', 'linewidth', 3, 'DisplayName','E obs');

plot(log10(result_t0.loss), 'DisplayName','Loss t0');
plot(log10(result_t0.err_E), ':', 'linewidth', 3, 'DisplayName','E t0');
legend()





%% Decomposition
% L_est = ;
tic
result = L_decomposition_hamiltonian_kossakowski(result_prony.L_est, sysInfo.p, trueInfo);
L_decomp_time = toc
result = L_decomposition_hamiltonian_kossakowski(result_obs.L_est, sysInfo.p, trueInfo);


%% (One stage) ALS Hamiltonian and Kossakowski 
result = ALS_hamiltonian_kossakowski(all_rho_prony_pair_data, sysInfo, trueInfo);
