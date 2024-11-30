function [all_pair_data, derivative_err] = get_all_data_pair_full_state(all_rho, obsInfo, sysInfo, all_rho_obs, all_rho_prony)
% get all data pair for learning Lindbladian

[n, ~, ~, M] = size(all_rho);

%% Prony fit data pair
all_pair_data.prony = zeros(n, n, 2, M*(obsInfo.obs_len));

for t = 1:obsInfo.obs_len
    all_pair_data.prony(:, :, 1, (t-1)*M+1:t*M) = all_rho_obs(:, :, t, :);
    for i = 1:n
        for j = 1:n
            for m = 1:M
                all_pair_data.prony(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(obsInfo.tgrid(t));
            end
        end
    end
end

%% Observation data pair
% use finite difference for the approximation of the derivative 
% Time steps number less 1
all_pair_data.obs = zeros(n, n, 2, M*(obsInfo.obs_len));
for t = 1:obsInfo.obs_len
    all_pair_data.obs(:, :, 1, (t-1)*M+1:t*M) = all_rho_obs(:, :, t, :);
    all_pair_data.obs(:, :, 2, (t-1)*M+1:t*M) = (all_rho_obs(:, :, t+1, :) - all_rho_obs(:, :, t, :))/obsInfo.dt;
end



%% best data pair
all_pair_data.best = zeros(n, n, 2, M*(obsInfo.obs_len));
for t = 1:obsInfo.obs_len
    all_pair_data.best(:, :, 1, (t-1)*M+1:t*M) = all_rho(:, :, (t-1)*obsInfo.obs_gap+1, :);
    all_pair_data.best(:, :, 2, (t-1)*M+1:t*M) = (all_rho(:, :, (t-1)*obsInfo.obs_gap+2, :) - all_rho(:, :, (t-1)*obsInfo.obs_gap+1, :))/sysInfo.dt;
end



%% Prony fit data pair
all_pair_data.prony_t0 = all_pair_data.prony(:, :, :, 1:M);



%%
% % Prony fit data pair augmented
% all_rho_prony_pair_data_aug = zeros(n, n, 2, M*(sysInfo.steps+1));
% 
% for t = 1:obsInfo.obs_len+1
%     for i = 1:n
%         for j = 1:n
%             for m = 1:M
%                 all_rho_prony_pair_data_aug(i, j, 1, (t-1)*M+m) = all_rho_prony{i, j, m}.h(sysInfo.tgrid(t));
%                 all_rho_prony_pair_data_aug(i, j, 2, (t-1)*M+m) = all_rho_prony{i, j, m}.dh(sysInfo.tgrid(t));
%             end
%         end
%     end
% end



%%
derivative_err.prony = norm(all_pair_data.prony - all_pair_data.best, 'fro');
derivative_err.obs   = norm(all_pair_data.obs - all_pair_data.best, 'fro');
end
