function prony_result = Prony_fit_rho(all_rho, sysInfo)
% Use prony method to fit the trajectory of each entry of rho

%% 
n = sysInfo.n;
M = sysInfo.M;
dt = sysInfo.dt;


%% prony fitting of the density trajectories
obs_gap = 15;
obs_len = 24;
prony_p = 12;
result_h = cell(n, n, M);

for m = 1:M
    rho = all_rho(:, :, :, m);
    for ind_1 = 1:n
        for ind_2 = 1:n
            Prony_I.prony_p = prony_p;
            Prony_I.prony_N = obs_len;
            Prony_I.obs_dx  = obs_gap*dt;
            Prony_I.obs_ind = (0:obs_len)*obs_gap+1;
            Prony_I.obs_h_grid = squeeze(rho(ind_1, ind_2, Prony_I.obs_ind));
            Prony_I.polycoef_method = 'MP';
            Prony_I.weight_method = 'ID';
            Prony_I.root_normalization = 0;
            Prony_I.lambda_augmentation = 0;
            Prony_I.drop_0 = 0;

            result_h{ind_1, ind_2, m} = prony_method(Prony_I);
        end
    end
end
%%
ind_1 = 1;
ind_2 = 4;
Prony_I.obs_h_grid = squeeze(rho(ind_1, ind_2, (0:obs_len)*obs_gap+1));

figure;subplot(131);hold on;
plot(sysInfo.tgrid((0:obs_len)*obs_gap+1), real(Prony_I.obs_h_grid), '.', 'MarkerSize', 15)
plot(sysInfo.tgrid,                        real(squeeze(rho(ind_1, ind_2, :))));
plot(sysInfo.tgrid,                        real(result_h{ind_1, ind_2}.h(sysInfo.tgrid)));
plot(sysInfo.tgrid,                        real(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))));

plot(sysInfo.tgrid((0:obs_len)*obs_gap+1), imag(Prony_I.obs_h_grid), '.', 'MarkerSize', 15)
plot(sysInfo.tgrid,                        imag(squeeze(rho(ind_1, ind_2, :))));
plot(sysInfo.tgrid,                        imag(result_h{ind_1, ind_2}.h(sysInfo.tgrid)));
plot(sysInfo.tgrid,                        imag(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))));



xlim([0, sysInfo.tgrid(obs_len*obs_gap+1)])
legend('real obs','real noisy true', 'real fit', 'real true', 'imag obs','imag noisy true', 'imag fit', 'imag true')
title("\rho ")

subplot(132);hold on;
plot(sysInfo.tgrid((0:obs_len)*obs_gap+1), gradient(real(Prony_I.obs_h_grid), Prony_I.obs_dx), '.-', 'MarkerSize', 15)
plot(sysInfo.tgrid,                        gradient(real(squeeze(rho(ind_1, ind_2, :))),  sysInfo.dt));
plot(sysInfo.tgrid,                        real(result_h{ind_1, ind_2}.dh(sysInfo.tgrid)));
plot(sysInfo.tgrid,                        gradient(real(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))), sysInfo.dt));

plot(sysInfo.tgrid((0:obs_len)*obs_gap+1), gradient(imag(Prony_I.obs_h_grid), Prony_I.obs_dx), '.-', 'MarkerSize', 15)
plot(sysInfo.tgrid,                        gradient(imag(squeeze(rho(ind_1, ind_2, :))),  sysInfo.dt));
plot(sysInfo.tgrid,                        imag(result_h{ind_1, ind_2}.dh(sysInfo.tgrid)));
plot(sysInfo.tgrid,                        gradient(imag(squeeze(true_Info_plot.all_rho_true(ind_1, ind_2, :, 1))), sysInfo.dt));

xlim([0, sysInfo.tgrid(obs_len*obs_gap+1)])
legend('real obs','real noisy true', 'real fit', 'real true', 'imag obs','imag noisy true', 'imag fit', 'imag true')
title("\rho '")
% set(gcf,'LineWidth',3)

subplot(133);hold on;
eigL = eig(true_Info_plot.L_true);
xL = max(max(abs(real(eigL))), 1);
yL = max(max(abs(imag(eigL))), 1);
plot(eigL, '.', 'markersize', 20);
plot(result_h{ind_1, ind_2}.lam, '.', 'markersize', 20);
xlim([-xL, xL]); ylim([-yL, yL])
xline(0)
grid on;
title('Eigenvalues of L')

%%

figure;hold on;
for ind_1 = 1:sysInfo.n
    for ind_2 = 1:sysInfo.n
        plot(result_h{ind_1, ind_2}.lam, '.', 'markersize', 5);
    end
end


eigL = eig(true_Info_plot.L_true);
xL = max(max(abs(real(eigL))), 1);
yL = max(max(abs(imag(eigL))), 1);
xL = 5;
yL = 5;
plot(eigL, '.', 'markersize', 20);
xlim([-xL, xL]); ylim([-yL, yL])
xline(0)
grid on;
title('Eigenvalues of L')

end



