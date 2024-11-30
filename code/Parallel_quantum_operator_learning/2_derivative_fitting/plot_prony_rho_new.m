function plot_prony_rho_new(rho_prony_func, rho, rho_obs, sysInfo, obsInfo)
%% Plot an example of prony fitting result

tgrid       = sysInfo.tgrid;
% rho         = squeeze(all_rho(ind_1, ind_2, :, m));
d_rho       = gradient(rho, sysInfo.dt);

obs_tgrid   = obsInfo.tgrid;
% rho_obs     = squeeze(all_rho_obs(ind_1, ind_2, :, m));
d_rho_obs   = gradient(rho_obs, obsInfo.dt);

rho_prony   = rho_prony_func.h(tgrid);
d_rho_prony = gradient(rho_prony, sysInfo.dt);



figure;subplot(121);hold on;
plot(tgrid, real(rho), 'LineWidth', 5);
plot(tgrid, real(rho_prony), 'LineWidth', 2);
plot(obs_tgrid, real(rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 1);
xline(obsInfo.T)
plot(tgrid, imag(rho), 'LineWidth', 5);
plot(tgrid, imag(rho_prony), 'LineWidth', 2);
plot(obs_tgrid, imag(rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 1);
xlim([0, obsInfo.T*1.5])
title('rho fitting')
legend('true', 'prony', 'obs')

subplot(122);hold on;
plot(tgrid, real(d_rho), 'LineWidth', 5);
plot(tgrid, real(d_rho_prony), 'LineWidth', 2);
plot(obs_tgrid, real(d_rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 1);
xline(obsInfo.T)
plot(tgrid, imag(d_rho), 'LineWidth', 5);
plot(tgrid, imag(d_rho_prony), 'LineWidth', 2);
plot(obs_tgrid, imag(d_rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 1);
xlim([0, obsInfo.T*1.5])
title('rho derivative fitting')
legend('true', 'prony', 'obs')

end
