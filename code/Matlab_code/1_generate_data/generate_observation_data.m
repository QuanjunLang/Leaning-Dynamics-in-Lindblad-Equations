function [all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo)
% generate observation data

obs_std = obsInfo.obs_std;
obs_gap = obsInfo.obs_gap;
obs_len = obsInfo.obs_len;

obsInfo.dt      = sysInfo.dt*obs_gap;
obsInfo.T       = obsInfo.dt*obs_len;
obsInfo.tgrid   = 0:obsInfo.dt:obsInfo.T;
obsInfo.obs_ind = (0:obs_len)*obs_gap+1;

%%
if sysInfo.FULL_STATE
    all_rho_obs = all_rho(:, :, obsInfo.obs_ind, :);
    % all_rho_obs = all_rho_obs + max(abs(all_rho_obs), [], 'all')*randn(size(all_rho_obs))*obs_std;
    all_rho_obs = all_rho_obs + mean(abs(all_rho_obs), 'all')*randn(size(all_rho_obs))*obs_std;
else
    all_rho_obs = all_rho(:, obsInfo.obs_ind, :);
    all_rho_obs = all_rho_obs + mean(abs(all_rho_obs), 'all')*randn(size(all_rho_obs))*obs_std;
end




end
