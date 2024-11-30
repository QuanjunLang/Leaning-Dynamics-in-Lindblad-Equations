clc
close all
clear all
rng(0)

addPaths

%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n       = 3;            % 
sysInfo.M       = 10;            % number of independent trajectories
sysInfo.dt      = 0.0001;        % true data generation time grid
sysInfo.p       = 2;            % number of jump operators
sysInfo.steps   = 10000;
sysInfo = update_sys(sysInfo);

sysInfo.PAPER_FIG_DIR = 'figure';
%%
[all_rho, trueInfo] = generate_data(sysInfo);

%% decomposition of the true L into H and K

% L_decomp_true = L_decomposition_hamiltonian_kossakowski(trueInfo.L_true, sysInfo.p, trueInfo);
%% generate observational data

obsInfo.obs_std = 1e-5;
obsInfo.obs_gap = 1000;
obsInfo.obs_len = 10;


[all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);


%% Demo
% plot_sample_traj(sysInfo)


%% Prony fitting test
all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);

ind_1 = 1;
ind_2 = 2;
m = 1;

plot_prony_rho(all_rho_prony, all_rho, all_rho_obs, sysInfo, obsInfo, ind_1, ind_2, m)
% plot_all_prony_modes(all_rho_prony, sysInfo, trueInfo)



%% Two stage, get L using Prony fitted data
[all_pair_data, derivative_err] = get_all_data_pair_full_state(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony);

result_obs = multi_mat_als(all_pair_data.obs, trueInfo.r_true, 'true_para', trueInfo);
result_best = multi_mat_als(all_pair_data.best, trueInfo.r_true, 'true_para', trueInfo);
result_prony = multi_mat_als(all_pair_data.prony, trueInfo.r_true, 'true_para', trueInfo);
result_prony_t0 = multi_mat_als(all_pair_data.prony_t0, trueInfo.r_true, 'true_para', trueInfo);
%%
figure;

hold on; grid on

plot(log10(result_obs.loss), 'DisplayName','Loss obs', 'color', color1);
plot(log10(result_obs.err_E), ':', 'linewidth', 3, 'DisplayName','E obs', 'color', color1);

plot(log10(result_best.loss), 'DisplayName','Loss best', 'color', color2);
plot(log10(result_best.err_E), ':', 'linewidth', 3, 'DisplayName','E best', 'color', color2);

plot(log10(result_prony.loss), 'DisplayName','Loss Prony', 'color', color3);
plot(log10(result_prony.err_E), ':', 'linewidth', 3, 'DisplayName','E Prony', 'color', color3);

plot(log10(result_prony_t0.loss), 'DisplayName','Loss Prony t_0', 'color', color4);
plot(log10(result_prony_t0.err_E), ':', 'linewidth', 3, 'DisplayName','E Prony t_0', 'color', color4);


legend()
fontsize(16,"points")
title('Estimation of L in Frobinius norm')

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

exportgraphics(gcf, 'n_6_L_loss.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);


%% Decomposition
result_obs.decomposition    = L_decomposition_hamiltonian_kossakowski(result_obs.L_est, sysInfo.p, trueInfo, 'ttl', 'obs', 'plotON', 0);
result_prony.decomposition  = L_decomposition_hamiltonian_kossakowski(result_prony.L_est, sysInfo.p, trueInfo, 'ttl', 'Prony', 'plotON', 0);
result_best.decomposition   = L_decomposition_hamiltonian_kossakowski(result_best.L_est, sysInfo.p, trueInfo, 'ttl', 'best', 'plotON', 0);
result_prony_t0.decomposition  = L_decomposition_hamiltonian_kossakowski(result_prony_t0.L_est, sysInfo.p, trueInfo, 'ttl', 'Prony', 'plotON', 0);



figure;

hold on; grid on

plot(log10(result_obs.decomposition.h_err), 'DisplayName','Hamiltonian error obs', 'color', color1);
plot(log10(result_obs.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error obs', 'color', color1);

plot(log10(result_best.decomposition.h_err), 'DisplayName','Hamiltonian error best', 'color', color2);
plot(log10(result_best.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error best', 'color', color2);


plot(log10(result_prony.decomposition.h_err), 'DisplayName','Hamiltonian error Prony', 'color', color3);
plot(log10(result_prony.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error Prony', 'color', color3);

plot(log10(result_prony_t0.decomposition.h_err), 'DisplayName','Hamiltonian error Prony t_0', 'color', color4);
plot(log10(result_prony_t0.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error Prony t_0', 'color', color4);

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


%% Trajactory prediction error
tgrid       = linspace(0, 50, 100);

trueInfo.pred_traj = generate_data_prediction(tgrid, trueInfo.H_true, trueInfo.C_true);
result_best.pred_traj = generate_data_prediction(tgrid, result_best.decomposition.H, result_best.decomposition.C);
result_obs.pred_traj = generate_data_prediction(tgrid, result_obs.decomposition.H, result_obs.decomposition.C);
result_prony.pred_traj = generate_data_prediction(tgrid, result_prony.decomposition.H, result_prony.decomposition.C);


[n, ~, TN, M] = size(trueInfo.pred_traj);

ind1 = randi(n);
ind2 = randi(n);
indm = randi(M);
figure;
subplot(121);hold on;
plot(tgrid, real(squeeze(trueInfo.pred_traj(ind1, ind2, :, indm))), 'Color',color1, 'DisplayName','true')
plot(tgrid, real(squeeze(result_best.pred_traj(ind1, ind2, :, indm))), 'Color',color2, 'DisplayName','best')
plot(tgrid, real(squeeze(result_obs.pred_traj(ind1, ind2, :, indm))), 'Color',color3, 'DisplayName','obs')
plot(tgrid, real(squeeze(result_prony.pred_traj(ind1, ind2, :, indm))), 'Color',color4, 'DisplayName','prony')

subplot(122);hold on;
plot(tgrid, imag(squeeze(trueInfo.pred_traj(ind1, ind2, :, indm))), 'Color',color1, 'DisplayName','true')
plot(tgrid, imag(squeeze(result_best.pred_traj(ind1, ind2, :, indm))), 'Color',color2, 'DisplayName','best')
plot(tgrid, imag(squeeze(result_obs.pred_traj(ind1, ind2, :, indm))), 'Color',color3, 'DisplayName','obs')
plot(tgrid, imag(squeeze(result_prony.pred_traj(ind1, ind2, :, indm))), 'Color',color4, 'DisplayName','prony')
legend()







%% (One stage) ALS Hamiltonian and Kossakowski 

decomp_result_one_stage_obs = ALS_hamiltonian_kossakowski(all_pair_data.obs, sysInfo, sysInfo.p, trueInfo);
decomp_result_one_stage_best = ALS_hamiltonian_kossakowski(all_pair_data.best, sysInfo, sysInfo.p, trueInfo);
decomp_result_one_stage_prony = ALS_hamiltonian_kossakowski(all_pair_data.prony, sysInfo, sysInfo.p, trueInfo);
decomp_result_one_stage_prony_t0 = ALS_hamiltonian_kossakowski(all_pair_data.prony_t0, sysInfo, sysInfo.p, trueInfo);


%%
lnwd = 1.5;
figure;hold on;grid on;
semilogy(decomp_result_one_stage_obs.h_err, 'DisplayName','Hamiltonian error obs', 'color', color1);
semilogy(decomp_result_one_stage_obs.K_err, ':', 'linewidth', 3, 'DisplayName','Kossa error obs', 'color', color1);

semilogy(decomp_result_one_stage_best.h_err, 'DisplayName','Hamiltonian error best', 'color', color2);
semilogy(decomp_result_one_stage_best.K_err, ':', 'linewidth', 3, 'DisplayName','Kossa error best', 'color', color2);

semilogy(decomp_result_one_stage_prony.h_err, 'DisplayName','Hamiltonian error prony', 'color', color3);
semilogy(decomp_result_one_stage_prony.K_err, ':', 'linewidth', 3, 'DisplayName','Kossa error prony', 'color', color3);

semilogy(decomp_result_one_stage_prony_t0.h_err, 'DisplayName','Hamiltonian error prony t_0', 'color', color4);
semilogy(decomp_result_one_stage_prony_t0.K_err, ':', 'linewidth', 3, 'DisplayName','Kossa error prony t_0', 'color', color4);



set(gca, 'YScale', 'log')
legend()
