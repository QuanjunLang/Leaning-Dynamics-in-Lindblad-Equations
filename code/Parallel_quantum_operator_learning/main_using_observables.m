clc
close all
clear all
rng(1)

addPaths
%% TODO
% Observable O and rho expectation equals to
% sum(O.' .* rho, 'all')
% which is the same as 
% vec(O')'*vec(rho)

%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n       = 4;            % dimension of rho
sysInfo.M       = 20;            % number of independent trajectories
sysInfo.dt      = 0.0001;       % true data generation time grid
sysInfo.p       = 4;            % number of jump operators
sysInfo.N_o     = 10;            % number of obserables
                                % set to n^2 indecate a full state observation                             
sysInfo.steps   = 10000;        % total number of time steps
sysInfo = update_sys(sysInfo);

[all_rho, trueInfo, observableInfo] = generate_data(sysInfo);

%% decomposition of the true L into H and K
% L_decomp_true = L_decomposition_hamiltonian_kossakowski(trueInfo.L_true, sysInfo.p, trueInfo);


%% generate observational data
obsInfo.obs_std = 0;
obsInfo.obs_gap = 1000;
obsInfo.obs_len = 10;

assert(obsInfo.obs_len*obsInfo.obs_gap <= sysInfo.steps, 'observation setting wrong')
[all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);


%% Prony fitting test
all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);
plot_prony_rho_one_sample(sysInfo, obsInfo, all_rho_prony, all_rho, all_rho_obs)


%% Construct input data for ALS with observables
[all_expect_derivative_0, derivative_err] = get_all_data_pair_observable(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony);


%% Two stage with observables

% result = multi_mat_als_observable(all_expect_derivative_0.best, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo);




%% Two stage
% get L using best data
result_best = multi_mat_als_observable(all_expect_derivative_0.best, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo, 'plotON', 0);

% get L using Prony fitted data
result_prony = multi_mat_als_observable(all_expect_derivative_0.prony, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo, 'plotON', 0);

% get L using finite difference derivative
result_obs = multi_mat_als_observable(all_expect_derivative_0.obs, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo, 'plotON', 0);



%%
figure;

hold on; grid on
plot(log10(result_prony.loss), 'DisplayName','Loss Prony', 'color', color1);
plot(log10(result_prony.err_E), ':', 'linewidth', 3, 'DisplayName','E Prony', 'color', color1);

plot(log10(result_obs.loss), 'DisplayName','Loss obs', 'color', color2);
plot(log10(result_obs.err_E), ':', 'linewidth', 3, 'DisplayName','E obs', 'color', color2);

% plot(log10(result_t0.loss), 'DisplayName','Loss t0', 'color', color3);
% plot(log10(result_t0.err_E), ':', 'linewidth', 3, 'DisplayName','E t0', 'color', color3);

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
% result_t0.decomposition     = L_decomposition_hamiltonian_kossakowski(result_t0.L_est, sysInfo.p, trueInfo, 'ttl', 't0', 'plotON', 0);
result_best.decomposition   = L_decomposition_hamiltonian_kossakowski(result_best.L_est, sysInfo.p, trueInfo, 'ttl', 'best', 'plotON', 0);

%%
figure;

hold on; grid on
plot(log10(result_prony.decomposition.h_err), 'DisplayName','Hamiltonian error Prony', 'color', color1);
plot(log10(result_prony.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error Prony', 'color', color1);

plot(log10(result_obs.decomposition.h_err), 'DisplayName','Hamiltonian error obs', 'color', color2);
plot(log10(result_obs.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error obs', 'color', color2);

% plot(log10(result_t0.decomposition.c_err), 'DisplayName','c error t0', 'color', color3);
% plot(log10(result_t0.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','K error t0', 'color', color3);

plot(log10(result_best.decomposition.h_err), 'DisplayName','Hamiltonian error best', 'color', color4);
plot(log10(result_best.decomposition.K_err), ':', 'linewidth', 3, 'DisplayName','Kossa error best', 'color', color4);

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




%%

% % K = result_best.decomposition.K;
% p = sysInfo.p;
% n = sysInfo.n;
% F = trueInfo.F;
% K = trueInfo.Kossakowski;
% C = cell(p, 1);
% [U, S, ~] = svd(K);
% US = U*sqrt(S);
% for i = 1:p
%     C{i} = zeros(n, n);
%     for j = 1:n^2-1
%         C{i} = C{i} + US(j, i)*F{j};
%     end
% end
% 
% C_true = trueInfo.C_true;
% 
% C_true{1}
% C{1}
% 
% 
% 
% C_true{1}./C{1}
% C_true{1}./result_best.decomposition.C{1}

% C{1}./result_best.decomposition.C{1}



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


