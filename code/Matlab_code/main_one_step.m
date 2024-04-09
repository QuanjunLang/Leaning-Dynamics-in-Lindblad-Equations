clc
close all
clear all
rng(1)
% pyenv(Version='/opt/anaconda3/envs/py3.11/bin/python')

%%
addPaths
% sysInfo.n       = 4;
% sysInfo.steps   = 20;
% sysInfo.dt      = 0.01;
% sysInfo.p       = 1;
% sysInfo.M       = 500;
% sysInfo = update_sys(sysInfo);


%%
sysInfo.n       = 4;
sysInfo.dt      = 0.01;
sysInfo.p       = 1;
sysInfo.M       = 12;
sysInfo.obs_std = 0.02;

sysInfo.tgrid   = [0, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10];
% Time grid for which we have observations
sysInfo.steps   = length(sysInfo.tgrid)-1;
sysInfo.L       = length(sysInfo.tgrid);
sysInfo.T       = 50;

%% Demo
% plot sample trajectories
plotON = 1;
if plotON


    sysInfo_plot        = sysInfo;
    sysInfo_plot.M      = 1;
    sysInfo_plot.dt     = 0.1;
    sysInfo_plot.steps  = sysInfo.T/sysInfo_plot.dt;
    sysInfo_plot = update_sys(sysInfo_plot);
    n = sysInfo_plot.n;
    tgrid = sysInfo_plot.tgrid;

    [sysInfo_plot, all_rho, true_Info_plot] = generate_data(sysInfo_plot, 'plotON', 1);
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




%% Prony fitting test
prony_fitting_test;

% prony_result = Prony_fit_rho(all_rho, sysInfo)






%% Learn the channel operators for different time difference
[sysInfo, all_rho, trueInfo] = generate_data(sysInfo);

all_rho_diff_time   = cell(sysInfo.steps, 1);
result_channel      = cell(sysInfo.steps, 1);

for t = 1:sysInfo.steps
    temp = zeros(sysInfo.n, sysInfo.n, 2, sysInfo.M);
    temp(:, :, 1, :) = all_rho(:, :, 1, :);
    % temp(:, :, 2, :) = all_rho(:, :, 1+t, :);
    temp(:, :, 2, :) = (all_rho(:, :, 1+t, :) - all_rho(:, :, 1, :))/sysInfo.tgrid(t+1);
    all_rho_diff_time{t} = temp;
end

for t = 1:sysInfo.steps
    % result_channel{t} = multi_mat_als_update_rank(all_rho_diff_time{t}, sysInfo.n^2-1, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 1);
    result_channel{t} = multi_mat_als(all_rho_diff_time{t}, sysInfo.n^2-11, 'true_para', trueInfo, 'rank_threshold', 1e-5, 'plotON', 0);
end

%%
figure;
for t = 1:8
    subplot(2, 4, t);hold on; grid on
    plot(log10(result_channel{t}.loss), 'DisplayName','Loss');
    plot(log10(result_channel{t}.err_E), ':', 'linewidth', 3, 'DisplayName','E');
    ttl = ['t = ', num2str(sysInfo.tgrid(t+1))];
    title(ttl)
    if t == 1
        legend()
    end
end
%% E is the low rank form
all_ES = zeros(sysInfo.n^2, sysInfo.steps);
all_EU = zeros(sysInfo.n^2, sysInfo.n^2, sysInfo.steps);
all_E  = zeros(sysInfo.n^2, sysInfo.n^2, sysInfo.steps);
for t = 1:sysInfo.steps
    % rank(result_channel{t}.E_est)
    % eig(result_channel{t}.S_est)
    all_E(:, :, t) = result_channel{t}.E_est;
    [U, S, V] = svd(all_E(:, :, t));
    all_ES(:, t) = diag(S);
    all_EU(:, :, t) = U;
end


%% L is the full rank form
all_LS = zeros(sysInfo.n^2, sysInfo.steps);
all_LU = zeros(sysInfo.n^2, sysInfo.n^2, sysInfo.steps);
all_L  = zeros(sysInfo.n^2, sysInfo.n^2, sysInfo.steps);
for t = 1:sysInfo.steps
    all_L(:, :, t) = result_channel{t}.L_est;
    [U, S, V] = svd(all_L(:, :, t));
    all_LS(:, t) = diag(S);
    all_LU(:, :, t) = U;
end


%% compare the estimation error of the operator E
err_E = zeros(sysInfo.steps, 1);
for t = 1:sysInfo.steps
    err_E(t) = norm(all_E(:, :, t) - trueInfo.E_true, 'fro');
end


%% compare the estimation error of the operator L
err_L = zeros(sysInfo.steps, 1);
for t = 1:sysInfo.steps
    err_L(t) = norm(all_L(:, :, t) - trueInfo.L_true, 'fro');
end


figure;
subplot(221);hold on;grid on;
plot(log10(sysInfo.tgrid(2:end)), log10(all_ES)')
xlabel('log_{10} time')
title('log_{10} E eigen values')


subplot(222);hold on;grid on;
plot(log10(sysInfo.tgrid(2:end)), log10(all_LS)')
xlabel('log_{10} time')
title('log_{10} L eigen values')


subplot(223);hold on;
plot(log10(sysInfo.tgrid(2:end)), log10(err_E), 'DisplayName','Err E')
plot(log10(sysInfo.tgrid(2:end)), log10(err_L), 'o', 'DisplayName','Err L')
legend('Location','best')
title('log 10 error of E and L with different time step size')
xlabel('log10 of time step')

subplot(224);hold on;
plot(log10(all_ES(:, 3)), '.-', 'MarkerSize',10);
plot(log10(all_LS(:, 3)), '.-', 'MarkerSize',10);
title('SVD of estimated E and L')


set(gcf,'Position',[100 100 1000 600])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%%
% power_err = zeros(sysInfo.steps, 1);
% for t = 1:sysInfo.steps
%     power_err(t) = norm(all_L(:, :, 1)^t - all_L(:, :, t), 'fro');
% end
% 
% figure;plot(log10(sysInfo.tgrid(2:end)), power_err)
% 

%%

