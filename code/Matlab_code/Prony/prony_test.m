clc
close all
clear all
rng(4)
rng(5)
rng(6)
% rng(8)
%%
obs_dx = 0.1;
N = 100;
obs_xgrid = (0:obs_dx:(N-1)*obs_dx)';
obs_std = 0.0001;
p = 50;
dx = 0.01;
xgrid = (0:dx:20)';



h_grid_obs = mlf(3/2, 1, -abs(obs_xgrid).^(3/2), 10) + randn(size(obs_xgrid))*obs_std;
h_grid = mlf(3/2, 1, -abs(xgrid).^(3/2), 10);
%%
% I_rkhs = polynomial_method(h_grid_obs, p, obs_dx, 'RKHS');
% I_id = polynomial_method(h_grid_obs, p, obs_dx, 'ID');
I_poly = polynomial_method(h_grid_obs, p, obs_dx, 'LS');
I_mtpc = matrix_pencil(h_grid_obs, p, obs_dx);
I_true = polynomial_method(h_grid, p, dx, 'LS');

I_rkhs = polynomial_method(h_grid_obs, p, obs_dx, 'RKHS_h0');
I_id = polynomial_method(h_grid_obs, p, obs_dx, 'ID_h0');
a = 1;
%%


figure;hold on;
plot(xgrid, h_grid, '-', 'LineWidth', 4)
plot(obs_xgrid, h_grid_obs, '.', 'MarkerSize', 10, 'LineWidth', 3)

plot(xgrid, real(I_poly.f(xgrid)), '-.', 'MarkerSize', 2, 'LineWidth', 2);
plot(xgrid, real(I_mtpc.f(xgrid)), '--', 'MarkerSize', 2, 'LineWidth', 2);
plot(xgrid, real(I_rkhs.f(xgrid)), '--', 'MarkerSize', 2, 'LineWidth', 2);
plot(xgrid, real(I_id.f(xgrid)),   '--', 'MarkerSize', 2, 'LineWidth', 2);
plot(xgrid, real(I_true.f(xgrid)), '--', 'MarkerSize', 2, 'LineWidth', 2);

xlim([0, obs_xgrid(end)+3])
ylim([-1.5, 1.5])
legend('True', 'obs', 'poly', 'mtpc', 'rkhs', 'id', 'best')


%%
figure;hold on;
plot(I_poly.a, 'o');
plot(I_rkhs.a);
plot(I_id.a);
plot(I_true.a);

legend('poly', 'rkhs', 'id', 'true')
%%
figure

xrange = max([abs(I_poly.r); abs(I_id.r); abs(I_true.r); abs(I_rkhs.r); abs(I_mtpc.r);1.1]);

subplot(231);
scatter(real(I_poly.r), imag(I_poly.r), 130, '.'); hold on;circle(0, 0, 1)
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('poly')

subplot(232);
scatter(real(I_id.r), imag(I_id.r), 130, '.'); hold on;circle(0, 0, 1)
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('id')

subplot(233);
scatter(real(I_rkhs.r), imag(I_rkhs.r), 130, '.'); hold on;circle(0, 0, 1)
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('rkhs')

subplot(234);
scatter(real(I_true.r), imag(I_true.r), 130, '.'); hold on;circle(0, 0, 1)
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('true')

subplot(235);
scatter(real(I_mtpc.r), imag(I_mtpc.r), 130, '.'); hold on;circle(0, 0, 1)
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('matrix pencil')
set(gcf, 'Position', [2000, 2000, 1200, 800]);


%% Cauchy bound
figure;hold on;
plot_Cauchy_bound(I_poly)
plot_Cauchy_bound(I_rkhs)
plot_Cauchy_bound(I_id)
plot_Cauchy_bound(I_true)
legend('poly', 'rkhs', 'id', 'true')

%% 
partial_fraction_decomposition


%%




figure;hold on;
plot(xgrid, gradient(gradient(real(I_poly.f(xgrid)), dx)), '-.', 'MarkerSize', 2, 'LineWidth', 2);
plot(xgrid, gradient(gradient(h_grid', dx)), '--', 'MarkerSize', 2, 'LineWidth', 2);
