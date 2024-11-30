function plot_prony_rho(all_rho_prony, all_rho, all_rho_obs, sysInfo, obsInfo, ind_1, ind_2, m)
%% Plot an example of prony fitting result



color1 = [     0    0.4470    0.7410];
color2 = [0.8500    0.3250    0.0980];
color3 = [0.9290    0.6940    0.1250];
color4 = [0.4940    0.1840    0.5560];


tgrid       = sysInfo.tgrid;
rho         = squeeze(all_rho(ind_1, ind_2, :, m));
d_rho       = gradient(rho, sysInfo.dt);

obs_tgrid   = obsInfo.tgrid;
rho_obs     = squeeze(all_rho_obs(ind_1, ind_2, :, m));
d_rho_obs   = gradient(rho_obs, obsInfo.dt);

rho_prony   = all_rho_prony{ind_1, ind_2, m}.h(tgrid);
d_rho_prony = gradient(rho_prony, sysInfo.dt);


true_lnwd = 6;
true_lnst_real = '-';
true_lnst_imag = ':';

prony_lnwd = 2;
prony_lnst = '^-';
prony_gap = 60;
prony_mksz = 8;

obs_mksz = 7;
obs_lnwd = 2;
obs_lnst = '+-';



figure;subplot(121);hold on;
p1 = plot(tgrid, real(rho), true_lnst_real,'LineWidth', true_lnwd, 'Color', color1, 'DisplayName', 'True trajectory (real)');
p2 = plot(tgrid(1:prony_gap:end), real(rho_prony(1:prony_gap:end)), prony_lnst, 'LineWidth', prony_lnwd, 'MarkerSize', prony_mksz, 'Color', color2, 'DisplayName', 'Prony fitted');
p3 = plot(obs_tgrid, real(rho_obs), obs_lnst, 'MarkerSize', obs_mksz, 'LineWidth', obs_lnwd, 'Color', color3, 'DisplayName', 'Observation');

p4 = plot(tgrid, imag(rho), true_lnst_imag, 'LineWidth', true_lnwd, 'DisplayName', 'True trajectory (imag)', 'Color', color1);
plot(tgrid(1:prony_gap:end), imag(rho_prony(1:prony_gap:end)), prony_lnst, 'LineWidth', prony_lnwd, 'MarkerSize', prony_mksz, 'DisplayName', 'Prony fitted', 'Color', color2);
plot(obs_tgrid, imag(rho_obs), obs_lnst, 'MarkerSize', obs_mksz, 'LineWidth', obs_lnwd, 'DisplayName', 'observed', 'Color', color3);
xlim([0, obsInfo.T])

title('Trajectory fitting')

xlabel('time t')


subplot(122);hold on;
plot(tgrid, real(d_rho), true_lnst_real, 'LineWidth', true_lnwd, 'Color', color1, 'DisplayName', 'true');
plot(tgrid(1:prony_gap:end), real(d_rho_prony(1:prony_gap:end)), prony_lnst, 'LineWidth', prony_lnwd, 'MarkerSize', prony_mksz, 'Color', color2, 'DisplayName', 'Prony fitted');
plot(obs_tgrid, real(d_rho_obs), obs_lnst, 'MarkerSize', obs_mksz, 'LineWidth', obs_lnwd, 'Color', color3, 'DisplayName', 'observed');


plot(tgrid, imag(d_rho), true_lnst_imag, 'LineWidth', true_lnwd, 'Color', color1, 'DisplayName', 'true');
plot(tgrid(1:prony_gap:end), imag(d_rho_prony(1:prony_gap:end)), prony_lnst, 'LineWidth', prony_lnwd, 'MarkerSize', prony_mksz, 'Color', color2, 'DisplayName', 'Prony fitted');
plot(obs_tgrid, imag(d_rho_obs), obs_lnst, 'MarkerSize', obs_mksz, 'LineWidth', obs_lnwd, 'Color', color3, 'DisplayName', 'observed');
xlim([0, obsInfo.T])


title('Trajectory derivative fitting')
xlabel('time t')


legend([p1, p4, p2, p3], 'Location','best')




%%%%%%%%%%%% save %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'Position',[100 100 1300 300])
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

figname = [sysInfo.PAPER_FIG_DIR, '/sample_trajectory_estimation.pdf'];
% saveas(gcf, figname);


exportgraphics(gcf, figname, 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution', 300);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

