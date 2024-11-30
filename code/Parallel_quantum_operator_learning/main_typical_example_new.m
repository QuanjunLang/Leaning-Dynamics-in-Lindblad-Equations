clc
close all
clear all
rng(0)

addPaths

%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n               = 16;            %
sysInfo.M               = 200;          % number of independent trajectories
sysInfo.dt              = 0.001;      % true data generation time grid
sysInfo.p               = 1;            % number of jump operators
sysInfo.steps           = 10;
sysInfo.channel_dt_rate = 10;



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
E_sub_blocks = trueInfo.E_sub_blocks;
L_sub_blocks = trueInfo.L_sub_blocks;
E_obs_blocks = trueInfo.E_obs_blocks;
L_obs_blocks = trueInfo.L_obs_blocks;

all_rho0    = observableInfo.rho0;



%%



%% single block test
ind_m = randi([1, sysInfo.M]);
ind_o = randi([1, sysInfo.N_o]);

rho0_O = all_rho(ind_o, 1, ind_m);
rho1_O = all_rho(ind_o, 2, ind_m);

O = observableInfo.O(:, :, ind_o);
rho0 = observableInfo.rho0(:, :, ind_m);

rho0_prime = (rho1_O - rho0_O)/sysInfo.dt;

L_error_vec = norm(vec(rho0)' * L' * vec(O) - rho0_prime);
E_error_vec = norm(vec(rho0)' * E' * vec(O) - rho1_O);

fprintf('\n\nError Using vectorization:\n')
fprintf('Lindbladian error, vectorization: %.18f \n', L_error_vec)
fprintf('Channel error, vectorization: %.18f \n', E_error_vec)


UV = kron(conj(O), rho0);
fprintf('\n\nError Using rearrangement:\n')
L_error_R = norm(sum(conj(RL).*(UV), 'all') - rho0_prime);
E_error_R = norm(sum(conj(RE).*(UV), 'all') - rho1_O);
fprintf('Lindbladian error, R: %.18f \n', L_error_R)
fprintf('Channel error, R: %.18f \n', E_error_R)


fprintf('\n\nBlock Error Using rearrangement:\n')
fprintf('n = %d, ind_o = %d\n', n, ind_o)
L_block_error_R = norm(sum(conj(L_obs_blocks{ind_o}).*rho0, 'all') - rho0_prime);
E_block_error_R = norm(sum(conj(E_obs_blocks{ind_o}).*rho0, 'all') - rho1_O);
fprintf('Lindbladian block error, R: %.18f \n', L_block_error_R)
fprintf('Channel block error, R: %.18f \n', E_block_error_R)


%%

ind_o = 33;

% b_E = squeeze(all_rho(ind_o, sysInfo.channel_dt_rate+1, :));
% [RE_est, outputInfo_E] = ALS(all_rho0, b_E, rank(E_obs_blocks{ind_o}), 'X_true', E_obs_blocks{ind_o}, 'debugON', 1);


b_L = squeeze(all_rho(ind_o, 2, :) - all_rho(ind_o, 1, :))/sysInfo.dt;
b_L_est = squeeze(sum(conj(L_obs_blocks{ind_o}).*all_rho0, [1,2]));
l = L_obs_blocks{ind_o};

[RL_est, outputInfo_L] = ALS(all_rho0, b_L_est, rank(L_obs_blocks{ind_o}), 'X_true', L_obs_blocks{ind_o}, 'debugON', 1);




%%
n = sysInfo.n;
M = sysInfo.M;

operator_A = zeros(n^2, M);
for m = 1:M
    operator_A(:, m) = vec(all_rho0(:, :, m)).';
end


%%
N_o = sysInfo.N_o;

all_L_blocks_svd = zeros(n, N_o);
all_L_blocks_rank = zeros(N_o, 1);
for ind_o = 1:N_o
    all_L_blocks_svd(:, ind_o) = svd(L_obs_blocks{ind_o});
    all_L_blocks_rank(ind_o) = rank(L_obs_blocks{ind_o});
end

all_E_blocks_svd = zeros(n, N_o);
all_E_blocks_rank = zeros(N_o, 1);
for ind_o = 1:N_o
    all_E_blocks_svd(:, ind_o) = svd(E_obs_blocks{ind_o});
    all_E_blocks_rank(ind_o) = rank(E_obs_blocks{ind_o});
end


figure;
subplot(121);hold on;
plot(all_L_blocks_rank, 'b')
plot(all_E_blocks_rank, 'r')
title('rank of subblocks of E and L')
legend('L', 'E')
xlabel('index of observables')
subplot(122);hold on;
plot(log10(all_L_blocks_svd), 'b')
plot(log10(all_E_blocks_svd), 'r')
title('log10 svd of subblocks of E and L')
xlabel('index of observables')