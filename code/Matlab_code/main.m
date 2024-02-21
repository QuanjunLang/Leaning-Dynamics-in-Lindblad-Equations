clc
close all
clear all
rng(1)
% pyenv(Version='/opt/anaconda3/envs/py3.11/bin/python')


sysInfo.n       = 2;
sysInfo.steps   = 1;
sysInfo.dt      = 0.00001;
sysInfo.p       = 4;
sysInfo.M       = 10;

sysInfo = update_sys(sysInfo);


%% Demo


% plot sample trajectories
plotON = 1;
if plotON
    figure;
    hold on;
    grid on;
    view(10, 10)

    sysInfo_plot        = sysInfo;
    sysInfo_plot.M      = 1;
    sysInfo_plot.dt     = 0.1;
    sysInfo_plot.steps  = 1000;
    sysInfo_plot = update_sys(sysInfo_plot);
    n = sysInfo_plot.n;
    tgrid = sysInfo_plot.tgrid;

    [sysInfo_plot, all_rho, ~] = generate_data(sysInfo_plot);
    rho = all_rho(:, :, :, 1);
    for i = 1:n
        for j = 1:n
            plot3(real(squeeze(rho(i, j, :))), imag(squeeze(rho(i, j, :))), tgrid, 'linewidth', 3, ...
                'DisplayName',['(', num2str(i), ',' num2str(j), ')']);

        end
    end
    legend()

end


%% Generate data
[sysInfo, all_rho, trueInfo] = generate_data(sysInfo);

C_true = trueInfo.C_true;
H_true = trueInfo.H_true;


%% Demo Learning Hermitian

H_est = learning_Hermitian(sysInfo, all_rho, C_true);
H_err = Hermitian_error(H_est, H_true);

fprintf('Hermitian estimation error is %.6f \n', H_err)
%% Demo Learnng Lindbladian
for ind = 1:sysInfo.p
    C_other = C_true;
    C_other(ind) = [];

    C_est = learning_Lindbladian_ind(sysInfo, all_rho, trueInfo.H_true, C_other, C_true, ind);
    C_err = norm(svd(C_true{ind}) - svd(C_est));
    fprintf('The No.%d Lindbladian estimation error is %.6f \n', ind, C_err)

end
%%
C_est = learning_Lindbladian_ind_one_step(sysInfo, all_rho, H_true, C_other, C_true{ind})
C_err = norm(svd(C_true{ind}) - svd(C_est))












