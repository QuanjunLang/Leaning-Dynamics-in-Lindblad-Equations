clc
close all
clear all
rng(0)

addPaths

try
    terminate(pe)
catch
    pause(0.01)
end

% copyright - Quanjun Lang, 2024
%% system settings
sysInfo.n               = 20;           %
sysInfo.M               = 180;          % number of independent trajectories
sysInfo.dt              = 0.000001;     % true data generation time grid
sysInfo.p               = 1;            % number of jumpoperators
sysInfo.steps           = 3;
sysInfo.channel_dt_rate = 3;
sysInfo.observable_option  = 'First_row_col_diag';


sysInfo.PAPER_FIG_DIR = 'figure';
%%
all_n       = 12:1:20;
all_M       = 80:20:250;


Num_M = length(all_M);
Num_n = length(all_n);
Num_samples = 3;

all_RL_est  = cell(Num_M, Num_n, Num_samples);
all_RL_Info = cell(Num_M, Num_n, Num_samples);

for i = 1:Num_M
    for j = 1:Num_n
        for k = 1:Num_samples
            sysInfo.M = all_M(i);
            sysInfo.n = all_n(j);
            sysInfo = update_sys(sysInfo);
        
            pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
            % terminate(pe)
        
            [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);
            bL = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
        
            [RL_est, RL_Info] = parallel_ALS_first_row(observableInfo.rho0, bL, trueInfo.rank_RL_true, 'X_true_sub_blocks', trueInfo.RL_sub_blocks, ...
                'X_true_obs_blocks', trueInfo.RL_obs_blocks, 'X_true', trueInfo.RL_true);
            fprintf('Parallel learning L, error = %.8f, time = %.2f seconds \n', RL_Info.rel_error_X, RL_Info.time)
        
            all_RL_est{i, j, k}   = RL_est;
            all_RL_Info{i, j, k}  = RL_Info;
        end
    end
end




%% collect data

all_error = zeros(Num_M, Num_n, Num_samples);
all_time  = zeros(Num_M, Num_n, Num_samples);

for i = 1:Num_M
    for j = 1:Num_n
        for k = 1:Num_samples
            all_error(i, j, k) = all_RL_Info{i, j, k}.rel_error_X;
            all_time(i, j, k)  = all_RL_Info{i, j, k}.time;
        end
    end
end

%%
figure;hold on
colors = colororder;
for i = 1:Num_M
    scatter(all_n, squeeze(log10(all_error(i, :, :))), 500, '.', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :))
end
% plot(all_n, all_n.*log(all_n)*trueInfo.rank_RL_true)
%%
figure;hold on
colors = colororder;
for i = 1:Num_n
    for k = 1:Num_samples
        if k==1
            vs = 'on';
        else
            vs = 'off';
        end
        lgd = ['n = ', num2str(all_n(i))];
        scatter(all_M, squeeze(log10(all_error(:, i, k))), 100, 'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), 'DisplayName', lgd,'HandleVisibility', vs)
    end
end


% Convert vector to string for legend
legendStrings = arrayfun(@num2str, all_n, 'UniformOutput', false);
% Add the legend using the vector
% legend(legendStrings, 'Location', 'best');
legend()

xlabel('M')
ylabel('log 10 relative error')


%%
figure;hold on;grid on;view(20, 20);
[XX, YY] = meshgrid(all_n, all_M);

for k = 1:Num_samples
    scatter3(XX, YY, squeeze(log10(all_error(:, :, k))), 500, 'k.');
end
% scatter3(XX, YY, all_error)


% Define the grid for x and z
[x, z] = meshgrid(-2*pi:0.1:2*pi, -2*pi:0.1:2*pi);

% Calculate y as a function of x
f = @(x) x.*log(x)*trueInfo.rank_RL_true;


y = f(XX);

ee = linspace(log10(min(all_error, [], 'all')), log10(max(all_error, [], 'all')), Num_M);
[~, ZZ] = meshgrid(all_n, ee);
surf(XX, y, ZZ)
ylabel('M')
xlabel('n')



%%
figure;hold on;grid on;view(0, 90);
all_percent = sum(all_error<0.01, 3)/Num_samples;
% plot(all_n, f(all_n), 'k', 'LineWidth',20)
% surf(XX, YY, all_percent)

% heatmap(all_percent, 'Colormap', jet)

imagesc(all_n, all_M, all_percent)
colormap('autumn')
plot(all_n, f(all_n), 'k', 'LineWidth', 5)
text(all_n(end-1), f(all_n(end-1))+20, 'M=rnlog(n)')

% s = 0.2;
% g = @(x) x.*log(x)*trueInfo.rank_RL_true./x.^s+30;
% plot(all_n, g(all_n), 'k', 'LineWidth', 5)
colorbar

xlabel('System size n')
ylabel('Trajectory number M')
title('Recovery rate')

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);

%% Time analysis
figure;hold on;grid on;view(20, 20);
[XX, YY] = meshgrid(all_n, all_M);

for k = 1:Num_samples
    scatter3(XX, YY, squeeze(log10(all_time(:, :, k))), 500, 'k.');
end

ylabel('M')
xlabel('n')



