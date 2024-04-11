function plot_prony_rho(all_rho_prony, all_rho, all_rho_obs, sysInfo, obsInfo, ind_1, ind_2, m)
%% Plot an example of prony fitting result


tgrid       = sysInfo.tgrid;
rho         = squeeze(all_rho(ind_1, ind_2, :, m));
d_rho       = gradient(rho, sysInfo.dt);

obs_tgrid   = obsInfo.tgrid;
rho_obs     = squeeze(all_rho_obs(ind_1, ind_2, :, m));
d_rho_obs   = gradient(rho_obs, obsInfo.dt);

rho_prony   = all_rho_prony{ind_1, ind_2, m}.h(tgrid);
d_rho_prony = gradient(rho_prony, sysInfo.dt);



figure;subplot(121);hold on;
plot(tgrid, real(rho), 'LineWidth', 2);
plot(tgrid, real(rho_prony), 'LineWidth', 2);
plot(obs_tgrid, real(rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 2);
xline(obsInfo.T)
plot(tgrid, imag(rho), 'LineWidth', 2);
plot(tgrid, imag(rho_prony), 'LineWidth', 2);
plot(obs_tgrid, imag(rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 2);

subplot(122);hold on;
plot(tgrid, real(d_rho), 'LineWidth', 2);
plot(tgrid, real(d_rho_prony), 'LineWidth', 2);
plot(obs_tgrid, real(d_rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 2);
xline(obsInfo.T)
plot(tgrid, imag(d_rho), 'LineWidth', 2);
plot(tgrid, imag(d_rho_prony), 'LineWidth', 2);
plot(obs_tgrid, imag(d_rho_obs), '.-', 'MarkerSize', 25, 'LineWidth', 2);

% subplot(133);hold on;
% eigL = eig(true_Info_plot.L_true);
% xL = max(max(abs(real(eigL))), 1);
% yL = max(max(abs(imag(eigL))), 1);
% plot(eigL, '.', 'markersize', 20);
% plot(all_rho_prony{ind_1, ind_2}.lam, '.', 'markersize', 20);
% xlim([-xL, xL]); ylim([-yL, yL])
% xline(0)
% grid on;
% title('Eigenvalues of L')

% 
% Prony_I.obs_h_grid = squeeze(rho(ind_1, ind_2, (0:obs_len)*obs_gap+1));
% 
% figure;subplot(131);hold on;
% plot(obs_tgrid, real(Prony_I.obs_h_grid), '.', 'MarkerSize', 15)
% plot(sysInfo.tgrid,                        real(squeeze(rho(ind_1, ind_2, :))));
% plot(sysInfo.tgrid,                        real(all_rho_prony{ind_1, ind_2}.h(sysInfo.tgrid)));
% plot(sysInfo.tgrid,                        real(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))));
% 
% plot(obs_tgrid, imag(Prony_I.obs_h_grid), '.', 'MarkerSize', 15)
% plot(sysInfo.tgrid,                        imag(squeeze(rho(ind_1, ind_2, :))));
% plot(sysInfo.tgrid,                        imag(all_rho_prony{ind_1, ind_2}.h(sysInfo.tgrid)));
% plot(sysInfo.tgrid,                        imag(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))));
% 
% 
% 
% xlim([0, sysInfo.tgrid(obs_len*obs_gap+1)])
% legend('real obs','real noisy true', 'real fit', 'real true', 'imag obs','imag noisy true', 'imag fit', 'imag true')
% title("\rho ")
% 
% subplot(132);hold on;
% plot(obs_tgrid, gradient(real(Prony_I.obs_h_grid), Prony_I.obs_dx), '.-', 'MarkerSize', 15)
% plot(sysInfo.tgrid,                        gradient(real(squeeze(rho(ind_1, ind_2, :))),  sysInfo.dt));
% plot(sysInfo.tgrid,                        real(all_rho_prony{ind_1, ind_2}.dh(sysInfo.tgrid)));
% plot(sysInfo.tgrid,                        gradient(real(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))), sysInfo.dt));
% 
% plot(obs_tgrid, gradient(imag(Prony_I.obs_h_grid), Prony_I.obs_dx), '.-', 'MarkerSize', 15)
% plot(sysInfo.tgrid,                        gradient(imag(squeeze(rho(ind_1, ind_2, :))),  sysInfo.dt));
% plot(sysInfo.tgrid,                        imag(all_rho_prony{ind_1, ind_2}.dh(sysInfo.tgrid)));
% plot(sysInfo.tgrid,                        gradient(imag(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))), sysInfo.dt));
% 
% xlim([0, sysInfo.tgrid(obs_len*obs_gap+1)])
% legend('real obs','real noisy true', 'real fit', 'real true', 'imag obs','imag noisy true', 'imag fit', 'imag true')
% title("\rho '")
% % set(gcf,'LineWidth',3)
% 
% subplot(133);hold on;
% eigL = eig(true_Info_plot.L_true);
% xL = max(max(abs(real(eigL))), 1);
% yL = max(max(abs(imag(eigL))), 1);
% plot(eigL, '.', 'markersize', 20);
% plot(all_rho_prony{ind_1, ind_2}.lam, '.', 'markersize', 20);
% xlim([-xL, xL]); ylim([-yL, yL])
% xline(0)
% grid on;
% title('Eigenvalues of L')

%%

