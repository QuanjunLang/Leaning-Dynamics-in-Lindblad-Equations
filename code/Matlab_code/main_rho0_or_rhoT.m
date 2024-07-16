clc
close all
clear all
rng(0)



addPaths

%%
color1 = [     0    0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];
color4 = [0.4940    0.1840    0.5560];



%% TODO
% 1. derivative error
%   other method of fitting trajectories
%   number of observables
% 2. number of observations
%   


%% system settings
% loadON = 0;
% if loadON && exist('n_4.mat', 'file')
%     load('n_4.mat');
% else
% pe = pyenv(Version='/usr/bin/python3', ExecutionMode = 'OutOfProcess');
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)
sysInfo.n       = 6;            %
sysInfo.M       = 8;            % number of independent trajectories
sysInfo.dt      = 0.0001;        % true data generation time grid
sysInfo.p       = 3;            % number of jump operators

sysInfo.steps   = 10000;
sysInfo = update_sys(sysInfo);
tic
[all_rho, trueInfo] = generate_data(sysInfo);
sysInfo.time = toc;
% end

%% decomposition of the true L into H and K
L_decomp_true = L_decomposition_hamiltonian_kossakowski(trueInfo.L_true, sysInfo.p, trueInfo);


%% generate observational data
obsInfo.obs_std = 0;
obsInfo.obs_gap = 1000;
obsInfo.obs_len = 10;


[all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);

%% Prony fitting test
all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);

ind_1 = 4;
ind_2 = 2;
m = 1;

plot_prony_rho(all_rho_prony, all_rho, all_rho_obs, sysInfo, obsInfo, ind_1, ind_2, m)
% plot_all_prony_modes(all_rho_prony, sysInfo, trueInfo)
%% Learn the channel operators using observed derivatives and prony-fitted derivatives

[all_rho_pair, derivative_err] = get_all_data_pair(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony);

%% Two stage
% get L using Prony fitted data
result_prony = multi_mat_als(all_rho_pair.prony, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);

% get L using finite difference derivative
result_obs = multi_mat_als(all_rho_pair.obs, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);

% get L using time 0 only
result_t0 = multi_mat_als(all_rho_pair.prony_t0, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);

% get L using best data
result_best = multi_mat_als(all_rho_pair.best, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);


%%
figure;

hold on; grid on
plot(log10(result_prony.loss), 'DisplayName','Loss Prony', 'color', color1);
plot(log10(result_prony.err_E), ':', 'linewidth', 3, 'DisplayName','E Prony', 'color', color1);

plot(log10(result_obs.loss), 'DisplayName','Loss obs', 'color', color2);
plot(log10(result_obs.err_E), ':', 'linewidth', 3, 'DisplayName','E obs', 'color', color2);

plot(log10(result_t0.loss), 'DisplayName','Loss t0', 'color', color3);
plot(log10(result_t0.err_E), ':', 'linewidth', 3, 'DisplayName','E t0', 'color', color3);

plot(log10(result_best.loss), 'DisplayName','Loss best', 'color', color4);
plot(log10(result_best.err_E), ':', 'linewidth', 3, 'DisplayName','E best', 'color', color4);

legend()
fontsize(16,"points")
title('Estimation of L in Frobinius norm')

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

exportgraphics(gcf, 'n_6_L_loss.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);
%% Decomposition
result_prony.decomposition  = L_decomposition_hamiltonian_kossakowski(result_prony.L_est, sysInfo.p, trueInfo, 'ttl', 'Prony', 'plotON', 0);
result_obs.decomposition    = L_decomposition_hamiltonian_kossakowski(result_obs.L_est, sysInfo.p, trueInfo, 'ttl', 'obs', 'plotON', 0);
result_t0.decomposition     = L_decomposition_hamiltonian_kossakowski(result_t0.L_est, sysInfo.p, trueInfo, 'ttl', 't0', 'plotON', 0);
result_best.decomposition   = L_decomposition_hamiltonian_kossakowski(result_best.L_est, sysInfo.p, trueInfo, 'ttl', 'best', 'plotON', 0);

%%
figure;

hold on; grid on
plot(log10(result_prony.decomposition.c_err), 'DisplayName','c error Prony', 'color', color1);
plot(log10(result_prony.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','K error Prony', 'color', color1);

plot(log10(result_obs.decomposition.c_err), 'DisplayName','c error obs', 'color', color2);
plot(log10(result_obs.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','K error obs', 'color', color2);

plot(log10(result_t0.decomposition.c_err), 'DisplayName','c error t0', 'color', color3);
plot(log10(result_t0.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','K error t0', 'color', color3);

plot(log10(result_best.decomposition.c_err), 'DisplayName','c error best', 'color', color4);
plot(log10(result_best.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','K error best', 'color', color4);

legend()
fontsize(16,"points")
title('Estimation of Hamiltonian and Kossakowski')

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

exportgraphics(gcf, 'n_6_K_c_error.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);
% %% (One stage) ALS Hamiltonian and Kossakowski 
% result = ALS_hamiltonian_kossakowski(all_rho_prony_pair_data, sysInfo, trueInfo);
