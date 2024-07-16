function [all_rho_pair, derivative_err] = get_all_data_pair(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony)

[n, ~, ~, M] = size(all_rho);


% Prony fit data pair
all_rho_pair.prony = zeros(n, n, 2, M*(obsInfo.obs_len+1));

for t = 1:obsInfo.obs_len
    all_rho_pair.prony(:, :, 1, (t-1)*M+1:t*M) = all_rho_obs(:, :, t, :);
    for i = 1:n
        for j = 1:n
            for m = 1:M
                all_rho_pair.prony(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(obsInfo.tgrid(t));
            end
        end
    end
end

% Observation data pair
all_rho_pair.obs = zeros(n, n, 2, M*(obsInfo.obs_len+1));
for t = 1:obsInfo.obs_len
    all_rho_pair.obs(:, :, 1, (t-1)*M+1:t*M) = all_rho_obs(:, :, t, :);
    all_rho_pair.obs(:, :, 2, (t-1)*M+1:t*M) = (all_rho_obs(:, :, t+1, :) - all_rho_obs(:, :, t, :))/obsInfo.dt;
end


% best data pair
all_rho_pair.best = zeros(n, n, 2, M*(obsInfo.obs_len+1));
for t = 1:obsInfo.obs_len
    t0 = (t-1)*obsInfo.obs_gap + 1;
    all_rho_pair.best(:, :, 1, (t-1)*M+1:t*M) = all_rho(:, :, t0, :);
    all_rho_pair.best(:, :, 2, (t-1)*M+1:t*M) = (all_rho(:, :, t0+1, :) - all_rho(:, :, t0, :))/sysInfo.dt;
end


% Prony fit data pair augmented
all_rho_pair.prony_aug = zeros(n, n, 2, M*(sysInfo.steps+1));

for t = 1:obsInfo.obs_len
    for i = 1:n
        for j = 1:n
            for m = 1:M
                all_rho_pair.prony_aug(i, j, 1, (t-1)*M+m) = all_rho_prony{i, j, m}.h(sysInfo.tgrid(t));
                all_rho_pair.prony_aug(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(sysInfo.tgrid(t));
            end
        end
    end
end

% Prony fit data time 0
all_rho_pair.prony_t0 = zeros(n, n, 2, M);

t = 1;
for i = 1:n
    for j = 1:n
        for m = 1:M
            all_rho_pair.prony_t0(i, j, 1, (t-1)*M+m) = all_rho_prony{i, j, m}.h(sysInfo.tgrid(t));
            all_rho_pair.prony_t0(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(sysInfo.tgrid(t));
        end
    end
end


derivative_err.prony = norm(all_rho_pair.prony - all_rho_pair.best, 'fro');
derivative_err.obs   = norm(all_rho_pair.obs - all_rho_pair.best, 'fro');


end