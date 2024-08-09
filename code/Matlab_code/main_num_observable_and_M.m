clc
close all
clear all
rng(0)
addPaths

%%
sysInfo.n       = 6;
sysInfo.p       = 2;            % number of jump operators
sysInfo.dt      = 0.0001;       % true data generation time grid
sysInfo.steps   = 10000;        % total number of time steps


obsInfo.obs_std = 0;
obsInfo.obs_gap = 1000;
obsInfo.obs_len = 6;


all_No = 20:3:35;
all_M  = 6:2:12;


%% Using observables
all_result_observable = cell(length(all_No), length(all_M));
for i = 1:length(all_No)
    for j = 1:length(all_M)

        pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
        terminate(pe)

        sysInfo.N_o = all_No(i);
        sysInfo.M = all_M(j);
        sysInfo = update_sys(sysInfo);

        fprintf('Using observables, N_o = %d, M = %d\n', sysInfo.N_o, sysInfo.M)

        [all_rho, trueInfo, observableInfo] = generate_data(sysInfo);
        [all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);
        all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);
        [all_expect_derivative_0, ~] = get_all_data_pair_observable(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony);

        result_observable.best = multi_mat_als_observable(all_expect_derivative_0.best, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo, 'plotON', 0);
        result_observable.prony = multi_mat_als_observable(all_expect_derivative_0.prony, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo, 'plotON', 0);
        result_observable.obs = multi_mat_als_observable(all_expect_derivative_0.obs, observableInfo, trueInfo.r_true, 'trueInfo', trueInfo, 'plotON', 0);

        all_result_observable{i, j} = result_observable;
    end
end


%% Using state tomography

all_result_full = cell(1, length(all_M));
for j = 1:length(all_M)

    pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
    terminate(pe)

    sysInfo.N_o = sysInfo.n^2; % Full state tomography
    sysInfo.M = all_M(j);
    sysInfo = update_sys(sysInfo);

    fprintf('Full state tomography, M = %d\n', sysInfo.M)

    [all_rho, trueInfo] = generate_data(sysInfo);
    [all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);
    all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);
    [all_rho_pair, derivative_err] = get_all_data_pair(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony);

    result.prony = multi_mat_als(all_rho_pair.prony, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);
    result.obs = multi_mat_als(all_rho_pair.obs, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);
    result.t0 = multi_mat_als(all_rho_pair.prony_t0, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);
    result.best = multi_mat_als(all_rho_pair.best, trueInfo.r_true, 'true_para', trueInfo, 'plotON', 0);

    all_result_full{j} = result;
end



%% collect data
err.observable.best     = zeros(length(all_No), length(all_M));
err.observable.prony    = zeros(length(all_No), length(all_M));
err.observable.obs      = zeros(length(all_No), length(all_M));

err.full.best      = zeros(length(all_No), length(all_M));
err.full.prony      = zeros(length(all_No), length(all_M));
err.full.obs      = zeros(length(all_No), length(all_M));
err.full.t0      = zeros(length(all_No), length(all_M));

for i = 1:length(all_No)
    for j = 1:length(all_M)
        err.observable.best(i, j) = all_result_observable{i,j}.best.err_E(end);
        err.observable.prony(i, j) = all_result_observable{i,j}.prony.err_E(end);
        err.observable.obs(i, j) = all_result_observable{i,j}.obs.err_E(end);
    end
end

for j = 1:length(all_M)
    err.observable.best(i+1, j) = all_result_full{j}.best.err_E(end);
    err.observable.prony(i+1, j) = all_result_full{j}.prony.err_E(end);
    err.observable.obs(i+1, j) = all_result_full{j}.obs.err_E(end);
end

for j = 1:length(all_M)
    err.full.best(j) = all_result_full{j}.best.err_E(end);
    err.full.prony(j) = all_result_full{j}.prony.err_E(end);
    err.full.obs(j) = all_result_full{j}.obs.err_E(end);
    err.full.t0(j) = all_result_full{j}.t0.err_E(end);
end


%%

all_No_ext = [all_No, sysInfo.n^2];
addPaths

%% Fix M, plot N_o
figure;hold on;
ind_o = 1;
plot(all_M, err.observable.best(ind_o, :), 'DisplayName', 'Best')
plot(all_M, err.observable.prony(ind_o, :), 'DisplayName', 'Prony')
plot(all_M, err.observable.obs(ind_o, :), 'DisplayName', 'Obs')
set(gca, 'YScale', 'log', 'XScale', 'log')
legend()


%% Fix N_o, plot M
figure;
for ind_m = 1:length(all_M)

    subplot(length(all_M), 1, ind_m);hold on;
    % ind_m = 1;
    plot(all_No_ext, err.observable.best(:, ind_m), 'o-', 'DisplayName', 'Best', 'Color',color1)
    plot(all_No_ext, err.observable.prony(:, ind_m), 'o-', 'DisplayName', 'Prony', 'Color',color2)
    plot(all_No_ext, err.observable.obs(:, ind_m), 'o-', 'DisplayName', 'Obs', 'Color',color3)
    scatter(sysInfo.n^2, err.full.t0(ind_m), '*', 'MarkerEdgeColor', color4, 'DisplayName', 'Prony using T0')

    xlabel('Number of Observables');
    ttl = ['M = ', num2str(all_M(ind_m))];
    title(ttl);
    set(gca, 'YScale', 'log')

    if ind_m == 1
        legend('Location','southwest');
    end
    grid on;
end

sgtitle(['Total dimension of density matrix: N = ', num2str(sysInfo.n)])


fontsize(16, "points")
% title('Estimation of Hamiltonian and Kossakowski')

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [0, 0, 10, 12])
set(gcf, 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

exportgraphics(gcf, 'conv_M_No_N6_0720_L6.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);


%%

save('num_O_M_conv_M_small_0720_short_time.mat', 'err', 'all_M', 'all_No', 'sysInfo')