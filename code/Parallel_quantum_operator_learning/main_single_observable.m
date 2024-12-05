clc
close all
clear all
rng(0)

addPaths

%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 12;            % 
sysInfo.M               = 1320;          % number of independent trajectories
sysInfo.dt              = 0.00001;      % true data generation time grid
sysInfo.p               = 1;            % number of jump operators
sysInfo.steps           = 100;
sysInfo.channel_dt_rate = 100;



% sysInfo.observable_option  = 'First_row_col_diag';
% sysInfo.observable_option  = 'Full_state';
% sysInfo.observable_option  = 'Multiple_random_observables';  sysInfo.N_o = 7;
sysInfo.observable_option  = 'Single_random_observable';

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


%% Single observable ALS
% Single_observable_ALS_test


%%
M = sysInfo.M;
N = sysInfo.n;
all_O = observableInfo.O;
all_rho0 = observableInfo.rho0;

RE = trueInfo.RE_true;
RL = trueInfo.RL_true;

A = cell(M, 1);
for m = 1:M
    O = all_O(:, :, m);
    rho0 = all_rho0(:, :, m);
    A{m} = kron(conj(O), rho0);
end

A_mat = zeros(N^2, N^2, M);
for m = 1:M
    A_mat(:, :, m) = A{m};
end



% fprintf('Estimating E: ')
% b_E = squeeze(all_rho(sysInfo.channel_dt_rate+1, :))';
% r_E = trueInfo.rank_RE_true;
% [RE_est, outputInfo_E] = ALS(A_mat, b_E, r_E, 'X_true', RE, 'debugON', 1, 'operator_name', 'E');
% fprintf('Time: %.3f\n', outputInfo_E.time)


fprintf('Estimating L: ')
b_L = squeeze(all_rho(2, :) - all_rho(1, :))'/sysInfo.dt;
r_L = trueInfo.rank_RL_true;
[RL_est, outputInfo_L] = ALS(A_mat, b_L, r_L, 'X_true', RL, 'debugON', 1, 'operator_name', 'L');
fprintf('Time: %.3f\n', outputInfo_L.time)