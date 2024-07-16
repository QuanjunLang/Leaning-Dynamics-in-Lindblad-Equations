function [all_expect_derivative_0, derivative_err] = get_all_data_pair_observable(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony)

[N_o, ~, M] = size(all_rho);


% Prony fit data pair
all_expect_derivative_0.prony = zeros(N_o, M);
for ind_o = 1:N_o
    for m = 1:M
        all_expect_derivative_0.prony(ind_o, m) = all_rho_prony{ind_o, m}.dh(0);
    end
end

% Observation data pair
all_expect_derivative_0.obs = zeros(N_o, M);
for ind_o = 1:N_o
    for m = 1:M
        all_expect_derivative_0.obs(ind_o, m) = (all_rho_obs(ind_o, 2, m) - all_rho_obs(ind_o, 1, m))/obsInfo.dt;
    end
end


% best data pair
all_expect_derivative_0.best = zeros(N_o, M);
for ind_o = 1:N_o
    for m = 1:M
        all_expect_derivative_0.best(ind_o, m) = (all_rho(ind_o, 2, m) - all_rho(ind_o, 1, m))/sysInfo.dt;
    end
end



derivative_err.prony = norm(all_expect_derivative_0.prony - all_expect_derivative_0.best, 'fro');
derivative_err.obs   = norm(all_expect_derivative_0.obs - all_expect_derivative_0.best, 'fro');
end