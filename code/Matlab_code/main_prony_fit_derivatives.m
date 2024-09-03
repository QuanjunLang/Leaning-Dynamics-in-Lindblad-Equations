clc
close all
clear all
rng(0)

addPaths

%% system settings
pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)

sysInfo.n       = 4;            % 
sysInfo.M       = 8;            % number of independent trajectories
sysInfo.dt      = 0.005;        % true data generation time grid
sysInfo.p       = 3;            % number of jump operators

sysInfo.steps   = 1000;
sysInfo = update_sys(sysInfo);

[all_rho, trueInfo] = generate_data(sysInfo);

%% decomposition of the true L into H and K

L_decomp_true = L_decomposition_hamiltonian_kossakowski(trueInfo.L_true, sysInfo.p, trueInfo);
%% generate observational data

obsInfo.obs_std = 0;
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

% %% DEBUG
% figure;hold on;
% plot(abs(all_rho_prony{i, j, m}.dh(obsInfo.tgrid)))
% obs_drho = squeeze((all_rho_obs(i, j, 2:end, m) - all_rho_obs(i, j, 1:end-1, m))/obsInfo.dt);
% plot(abs(obs_drho), '-o');
% % plot(obs_drho);


%% Two stage, get L using Prony fitted data
tic
result_prony = multi_mat_als(all_rho_prony_pair_data, trueInfo.r_true, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);
L_time = toc;
% result_prony_aug = multi_mat_als(all_rho_prony_pair_data_aug, trueInfo.r_true, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);

%% Two stage, get L using finite difference derivative
result_obs = multi_mat_als(all_rho_obs_pair_data, trueInfo.r_true, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);
%%
figure;

% subplot(121);
hold on; grid on
plot(log10(result_prony.loss), 'DisplayName','Loss Prony');
plot(log10(result_prony.err_E), ':', 'linewidth', 3, 'DisplayName','E Prony');
% title('Prony fit derivatives')

% subplot(122);hold on; grid on
plot(log10(result_obs.loss), 'DisplayName','Loss obs');
plot(log10(result_obs.err_E), ':', 'linewidth', 3, 'DisplayName','E obs');
% ttl = ['obs derivatives'];
% title(ttl)

legend()
% subplot(133);hold on; grid on
% plot(log10(result_prony_aug.loss), 'DisplayName','Loss');
% plot(log10(result_prony_aug.err_E), ':', 'linewidth', 3, 'DisplayName','E');
% ttl = ['prony augmented derivatives'];
% title(ttl)





%% Decomposition
% L_est = ;
tic
result = L_decomposition_hamiltonian_kossakowski(result_prony.L_est, sysInfo.p, trueInfo);
L_decomp_time = toc
result = L_decomposition_hamiltonian_kossakowski(result_obs.L_est, sysInfo.p, trueInfo);


%% (One stage) ALS Hamiltonian and Kossakowski 
tic
result = ALS_hamiltonian_kossakowski(all_rho_prony_pair_data, sysInfo, trueInfo);
One_step_time = toc;

% 
% 
% 
% %% Decomposition of Liouvillian into Hamiltonian and Jump
% E = trueInfo.E_true;
% L = trueInfo.L_true;
% phi = trueInfo.phi_true;
% kappa = trueInfo.kappa_true;
% kappa_part = trueInfo.kappa_part_true;
% 
% K = rearrangement_R(kappa_part);
% P = rearrangement_R(phi);
% 
% norm(P - K - E)
% 
% %% 
% H = trueInfo.H_true;
% H_part = trueInfo.H_part_true;
% jump_part = trueInfo.jump_part_true;
% 
% 
% F = cell(n^2-1, 1);
% for l = 1:n-1
%     temp = zeros(n, 1);
%     temp(1:l) = 1;
%     temp(l+1) = -l;
%     F{l} = 1i/sqrt(l*(l+1))*(diag(temp));
% end
% l = l+1;
% for k = 1:n
%     for j = 1:k-1
%         e_jk = zeros(n, n);
%         e_jk(j, k) = 1;
% 
%         F{l} = (e_jk + e_jk')/sqrt(2);l = l+1;
%         F{l} = -1i*(e_jk - e_jk')/sqrt(2);l = l+1;
% 
%     end
% end
% 
% 
% 
% 
% E = zeros(n^2- 1, n^2-1);
% for i = 1:n^2-1
%     for j = 1:n^2-1
%         E(i, j) = trace(F{i}*F{j}');
%     end
% end
% 
% 
% 
% 
% G = zeros(n^2, n^2, n^2-1, n^2-1);
% 
% for i = 1:n^2-1
%     for j = 1:n^2-1
%         G(:, :, i, j) = kron(conj(F{i}), F{j}) - 0.5*kron((F{i}'*F{j}).', eye(n)) - 0.5*kron(eye(n), F{i}'*F{j});
%     end
% end
% 
% % AA = reshape(G, [n^2, n^2, (n^2-1)^2]);
% A = zeros(n^4, (n^2-1)^2);
% b = zeros(n^4, 1);
% 
% k = 1;
% for i = 1:n^2
%     for j = 1:n^2
%         A(k, :) = reshape(G(i, j, :, :), [], 1);
%         b(k) = jump_part(i, j);
%         k = k + 1;
%     end
% end
% 
% c = A\b;
% norm(A*c - b)
% C = reshape(c, [n^2-1, n^2-1]);
% 
% C = (C + C')/2;
% % C is the Kossakoski matrix 
% %%
% [U, S] = svd(C);
% U = U';
% u = U(1, :)*sqrt(S(1, 1));
% 
% C_true = trueInfo.C_true{1};
% C_est = zeros(size(C_true));
% 
% for i = 1:n^2 - 1
%     C_est = C_est + u(i)*F{i};
% end
% 
% % C_est
% % C_true
% 
% r = mean(C_true./C_est, 'all');
% abs(C_true./C_est)
% 
% % C_est_v = zeros(size(C_true));
% % v = u*r;
% % for i = 1:n^2 - 1
% %     C_est_v = C_est_v + v(i)*F{i};
% % end
% % C_true./C_est_v
% 
% %%
% [U, S] = svd(C);
% % U = U';
% u = U*sqrt(S);
% v = zeros(size(u));
% 
% C_true = trueInfo.C_true;
% for q = 1:length(C_true)
%     for i = 1:n^2 - 1
%         v(i, q) = trace(C_true{q}'*F{i});
%     end
% end
% 
% 
% 
% 
