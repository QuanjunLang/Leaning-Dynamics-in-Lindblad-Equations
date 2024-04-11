function plot_sample_traj(sysInfo)


sysInfo_plot        = sysInfo;
sysInfo_plot.M      = 1;
sysInfo_plot.dt     = 0.1;
sysInfo_plot.steps  = 300;
sysInfo_plot = update_sys(sysInfo_plot);
n = sysInfo_plot.n;
tgrid = sysInfo_plot.tgrid;

[all_rho, ~] = generate_data(sysInfo_plot, 'plotON', 1);
rho = all_rho(:, :, :, 1);

figure;
hold on;
grid on;
view(10, 10)
for i = 1:n
    for j = 1:n
        plot3(real(squeeze(rho(i, j, :))), imag(squeeze(rho(i, j, :))), tgrid, 'linewidth', 3, ...
            'DisplayName',['(', num2str(i), ',' num2str(j), ')']);
    end
end
% legend()

end


