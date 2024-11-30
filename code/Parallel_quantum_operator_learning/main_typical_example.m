clc
close all
clear all
rng(0)

addPaths

%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 3;            % 
sysInfo.M               = 80;          % number of independent trajectories
sysInfo.dt              = 0.0001;      % true data generation time grid
sysInfo.p               = 2;            % number of jump operators
sysInfo.steps           = 30;
sysInfo.channel_dt_rate = 30;



sysInfo.observable_option  = 'First_row_col_diag';
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
plot(log10(s'))
xline(trueInfo.r_true)
xline(rank(trueInfo.RE_true))
for i = 1:length(all_dt)
    plot(r(i), log10(s(i, r(i))), 'o');
end

legendStrings = arrayfun(@num2str, all_dt, 'UniformOutput', false);
legend([legendStrings, 'rank(L)'], 'Location', 'best');

xlabel('index of singular values')
ylabel('log_{10} of singular values')
title('The change of channel operator rank with different time scale')

%% generate observational data
% 
% obsInfo.obs_std = 1e-5;
% obsInfo.obs_gap = 1000;
% obsInfo.obs_len = 10;
% 
% [all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);

%% Prony fitting test
% all_rho_prony = Prony_fit_rho(all_rho_obs, obsInfo);
% plot_prony_rho_one_sample(sysInfo, obsInfo, all_rho_prony, all_rho, all_rho_obs)

%% Check the reshaping is correct
% test_reshaping_vectorization_check


%% Single observable ALS
% Single_observable_ALS_test

%% Multiple block ALS

L = trueInfo.L_true;
E = trueInfo.E_true;
RL = trueInfo.RL_true;
RE = trueInfo.RE_true;


size(observableInfo.O)
size(observableInfo.rho0)
size(all_rho)


%%

switch sysInfo.observable_option
    case 'Multiple_random_observables'
        ind_m = randi([1, sysInfo.M]);
        ind_o = randi([1, sysInfo.N_o]);
        rho0_O_temp = all_rho(ind_o, 1, ind_m);
        rho1_O_temp = all_rho(ind_o, 2, ind_m);

        O_temp = observableInfo.O(:, :, ind_o, ind_m);
        rho0_temp = observableInfo.rho0(:, :, ind_m);

    case 'Single_random_observable'
        ind_m = randi([1, sysInfo.M]);
        rho0_O_temp = all_rho(1, ind_m);
        rho1_O_temp = all_rho(2, ind_m);

        O_temp = observableInfo.O(:, :, ind_m);
        rho0_temp = observableInfo.rho0(:, :, ind_m);

    case 'First_row_col_diag'
        ind_m = randi([1, sysInfo.M]);
        ind_o = randi([1, sysInfo.N_o]);
        % ind_o = 1;
        rho0_O_temp = all_rho(ind_o, 1, ind_m);
        rho1_O_temp = all_rho(ind_o, 2, ind_m);

        O_temp = observableInfo.O(:, :, ind_o);
        rho0_temp = observableInfo.rho0(:, :, ind_m);

        
end


% E_blocks = cell(sysInfo.n, sysInfo.n);
% for i = 1:sysInfo.n
%     for j = 1:sysInfo.n
%         E_blocks{i, j} = trueInfo.RE_true(i:sysInfo.n:end, j:sysInfo.n:end);
%     end
% end



%% Old rearrangemnt

% E_blocks = cell(sysInfo.N_o, 1);
% for i = 1:sysInfo.n
%     E_blocks{i} = trueInfo.RE_true(i:sysInfo.n:end, i:sysInfo.n:end);
% end
% 
% for k = 1:sysInfo.n-1
%     E_blocks{i+k} = trueInfo.RE_true(1:sysInfo.n:end, k:sysInfo.n:end) + trueInfo.RE_true(k:sysInfo.n:end, 1:sysInfo.n:end);
% end
% 
% for l = 1:sysInfo.n-1
%     E_blocks{i+k+l} = -1i*trueInfo.RE_true(1:sysInfo.n:end, l:sysInfo.n:end) +1i*trueInfo.RE_true(l:sysInfo.n:end, 1:sysInfo.n:end);
% end





rho0_prime_temp = (rho1_O_temp - rho0_O_temp)/sysInfo.dt;


L_error_vec = norm(vec(rho0_temp)' * L' * vec(O_temp) - rho0_prime_temp);
E_error_vec = norm(vec(rho0_temp)' * E' * vec(O_temp) - rho1_O_temp);

fprintf('Lindbladian error, vectorization: %.18f \n', L_error_vec)
% fprintf('Lindbladian error, vectorization prony: %.18f \n', L_error_vec_prony)
fprintf('Channel error, vectorization: %.18f \n', E_error_vec)


UV = kron(conj(rho0_temp), O_temp);
L_error_R = norm(sum(conj(RL).*(UV), 'all') - rho0_prime_temp);
E_error_R = norm(sum(conj(RE).*(UV), 'all') - rho1_O_temp);
E_error_block = norm(sum(conj(E_blocks{ind_o}).*conj(rho0_temp), 'all') - rho1_O_temp)


fprintf('Lindbladian error, R: %.18f \n', L_error_R)
fprintf('Channel error, R: %.18f \n', E_error_R)



%%
all_rho0    = observableInfo.rho0;
all_O       = observableInfo.O;
[N, ~, M]   = size(all_rho0);
A = cell(M, 1);
for m = 1:M
    A{m} = conj(all_rho0(:, :, m));
end
ind_o = randi([1, sysInfo.N_o]);

b_E = squeeze(all_rho(ind_o, sysInfo.channel_dt_rate+1, :));
[RE_est, outputInfo_E] = ALS(A, b_E, trueInfo.r_true);


all_error_E = squeeze(sum((abs(outputInfo_E.all_M - conj(E_blocks{ind_o}))).^2, [1, 2]));


disp(ind_o)

figure;
plot(log10(all_error_E))