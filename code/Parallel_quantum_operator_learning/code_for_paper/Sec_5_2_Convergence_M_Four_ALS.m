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
sysInfo.n               = 8;           %
sysInfo.M               = 100;          % number of independent trajectories
sysInfo.dt              = 0.000001;     % true data generation time grid
sysInfo.p               = 3;            % number of jumpoperators
sysInfo.steps           = 1;
sysInfo.channel_dt_rate = 1;
sysInfo.observable_option  = 'Channel_first_row_col';

sysInfo.PAPER_FIG_DIR = 'figure';


pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)


%%
all_M       = 15:2:80;
% all_M = 80;

Num_M = length(all_M);
Num_samples = 20;


all_ALS_N2 = cell(Num_M, Num_samples);
all_ALS_P = cell(Num_M, Num_samples);
all_ALS_N = cell(Num_M, Num_samples);
all_ALS_I = cell(Num_M, Num_samples);


for i = 1:Num_M
    for k = 1:Num_samples
        
        %% Load parameters
        sysInfo.M = all_M(i);
        sysInfo = update_sys(sysInfo);
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('The %d-th M = %d out of %d total choice, %d-th sample out of %d\n', i, sysInfo.M, Num_M, k, Num_samples)
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        %% Generate data
        if contains(sysInfo.observable_option, 'Channel')
            % Channel learning
            [all_rho, trueInfo, observableInfo] = generate_data_channel(sysInfo, 'plotON', 1);

            b_mat = vec(all_rho);
            b_all = all_rho;

        else
            % Lindbladian learning
            [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);

            b_mat = vec(squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt);
            b_all = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;

        end



        %% Multiple block ALS
        n = sysInfo.n;
        M = sysInfo.M;
        N_o = sysInfo.N_o;

        K = trueInfo.E_true;
        RK = trueInfo.RE_true;
        rank_RK = trueInfo.rank_RE_true;
        RK_sub_blocks = trueInfo.RE_sub_blocks;
        all_rho0    = observableInfo.rho0;


        A_mat = zeros(n^2, n^2, M*N_o);
        l = 1;
        for m = 1:M
            rho0 = all_rho0(:, :, m);
            for j = 1:N_o
                O = observableInfo.O(:, :, j);
                A_mat(:, :, l) = kron(conj(O), rho0);
                l = l+1;
            end
        end


        %% settings 
        compute_dnorm = 0;

        %% Total ALS

        fprintf('\nJoint ALS for the entire matrix...\n')
        [RK_est_total, RK_Info_total] = ALS(A_mat, b_mat, rank_RK, 'X_true', RK, ...
            'Nesterov_beta', 0.5, 'loss_tol', 1e-6, 'rel_err_tol', 1e-7, 'debugON', 0, 'maxIter', 500);

        fprintf('Finished in %.2f seconds \n', RK_Info_total.time)

        RK_Info_total = error_analysis(RK_Info_total, RK_est_total, RK, 'displayON', 1, 'compute_dnorm', compute_dnorm);

        %% Learning in parallel for all first row blocks
        [RK_est, RK_Info] = parallel_ALS_first_row_new(all_rho0, all_rho, rank_RK, 'X_true_sub_blocks', RK_sub_blocks, 'X_true', RK, 'debugON', 0, 'displayON', 1, 'compute_dnorm', compute_dnorm, 'num_retry', 3, ...
            'Nesterov_beta', 0.5, 'ALS_loss_tol', 1e-9, 'ALS_rel_err_tol', 1e-9, 'ALS_debugON', 0);


        %% Learning the first row jointly

        [RK_est_first_row, RK_est_first_row_Info] = ALS_first_row_joint_learning(all_rho0, all_rho, rank_RK, 'X_true_sub_blocks', RK_sub_blocks, 'X_true', RK, 'num_retry', 3, 'debugON', 0, 'displayON', 1, 'compute_dnorm', compute_dnorm, ...
            'Nesterov_beta', 0.5, 'ALS_loss_tol', 1e-8, 'ALS_rel_err_tol', 1e-8, 'ALS_debugON', 0);


        %% Learning a random subset of the first row jointly
        [RK_est_first_row_subset, RK_est_first_row_Info_subset] = ALS_first_row_joint_learning_random_subset(all_rho0, all_rho, rank_RK, 'sub_ind_ratio', 0.2, 'X_true_sub_blocks', RK_sub_blocks, 'X_true', RK, 'num_retry', 3, 'debugON', 0, 'plotON', 0, 'displayON', 1, 'compute_dnorm', compute_dnorm, ...
            'Nesterov_beta', 0.5, 'ALS_loss_tol', 1e-8, 'ALS_rel_err_tol', 1e-8, 'ALS_debugON', 0);


        %% store the numbers
        all_ALS_N2{i, k} = RK_Info_total;
        all_ALS_P{i, k} = RK_Info;
        all_ALS_N{i, k} = RK_est_first_row_Info;
        all_ALS_I{i, k} = RK_est_first_row_Info_subset;


    end
end




%% collect data
Num_items = 2;

all_error = zeros(4, Num_items, Num_M, Num_samples);

for i = 1:Num_M
    for k = 1:Num_samples
    
        all_error(1, 1, i, k) = all_ALS_N2{i, k}.rel_error_fro_X;
        % all_error(1, 2, i, k) = all_ALS_N2{i, k}.error_diamond;
        all_error(1, 2, i, k) = all_ALS_N2{i, k}.time;
    
        try
            all_error(2, 1, i, k) = all_ALS_P{i, k}.rel_error_fro_X;
        catch 
            all_error(2, 1, i, k) = 100;
        end

        % all_error(2, 2, i, k) = all_ALS_P{i, k}.error_diamond;
        all_error(2, 2, i, k) = all_ALS_P{i, k}.time;

        try
            all_error(3, 1, i, k) = all_ALS_N{i, k}.rel_error_fro_X;
        catch
            all_error(3, 1, i, k) = 100;
        end
        % all_error(3, 2, i, k) = all_ALS_N{i, k}.error_diamond;
        all_error(3, 2, i, k) = all_ALS_N{i, k}.time;

        try
            all_error(4, 1, i, k) = all_ALS_I{i, k}.rel_error_fro_X;
        catch
            all_error(4, 1, i, k) = 100;
        end
        % all_error(4, 2, i, k) = all_ALS_I{i, k}.error_diamond;
        all_error(4, 2, i, k) = all_ALS_I{i, k}.time;
    end
end





%% successful recovery rate 
figure;hold on;grid on;
colors(1, :) = color1;
colors(2, :) = color2;
colors(3, :) = color3;
colors(4, :) = color4;

marker{1} = '+';
marker{2} = 'o';
marker{3} = '*';
marker{4} = 's';

lineWidth = 1;

lgd{1} = 'ALS-N^2';
lgd{2} = 'ALS-P';
lgd{3} = 'ALS-N';
lgd{4} = 'ALS-I';

for i = 1:4
    success_points = sum(squeeze(all_error(i, 1, :, :))<1e-5, 2);
    scatter(all_M, success_points, 48, colors(i, :), 'Marker', marker{i}, 'LineWidth', lineWidth, 'displayname', lgd{i});
end

legend('Location','best')
ylb = ['# error \leq 1e-05 /', num2str(Num_samples), ' trials'];
ylabel(ylb);
xlabel('Number of independent trials')




set(gcf, 'Position',  [100, 100, 700, 300]);
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
tightfig(gcf)
saveas(gcf, 'Sec_5_2_Four_method_rate.pdf')
saveas(gcf, '/Users/quanjunlang/Documents/GitHub/Learning-operators-in-Lindblad-Master-Equation-Paper/figures/Sec_5_2_Four_method_rate.pdf')
%% computation time
figure;hold on;
colors(1, :) = color1;
colors(2, :) = color2;
colors(3, :) = color3;
colors(4, :) = color4;

marker{1} = '+';
marker{2} = 'o';
marker{3} = '*';
marker{4} = 's';

lineWidth = 2;

lgd{1} = 'ALS-N^2';
lgd{2} = 'ALS-P';
lgd{3} = 'ALS-N';
lgd{4} = 'ALS-I';


for i = 1:4
    success_points = mean(squeeze(all_error(i, 2, :, :)), 2);
    % if i == 1
    %     yaxix left
    % else
    %     yaxis right
    % end
    scatter(all_M, success_points, 52, 'Color', colors(i, :), 'Marker', marker{i}, 'LineWidth', lineWidth, 'displayname', lgd{i});
end
set(gca, 'YScale', 'log')

legend('Location','best')
ylb = 'Median computation time (s)';
ylabel(ylb);
xlabel('Number of independent trials')

set(gcf, 'Position',  [100, 100, 700, 300]);
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
tightfig(gcf)
saveas(gcf, 'Sec_5_2_Four_method_time_1.pdf')
saveas(gcf, '/Users/quanjunlang/Documents/GitHub/Learning-operators-in-Lindblad-Master-Equation-Paper/figures/Sec_5_2_Four_method_time_1.pdf')
%% computation time (two figures)

colors(1, :) = color1;
colors(2, :) = color2;
colors(3, :) = color3;
colors(4, :) = color4;

marker{1} = '+';
marker{2} = 'o';
marker{3} = '*';
marker{4} = 's';

lineWidth = 2;

lgd{1} = 'ALS-N^2';
lgd{2} = 'ALS-P';
lgd{3} = 'ALS-N';
lgd{4} = 'ALS-I';


hfig = figure;
% subplot(211);

tiledlayout(2, 1, 'padding', 'compact', "TileSpacing","compact");
nexttile;
hold on;
% grid on;
i = 1;
success_points = median(squeeze(all_error(i, 2, :, :)), 2);
scatter(all_M, success_points, 48, colors(i, :), 'Marker', marker{i}, 'LineWidth', lineWidth, 'displayname', lgd{i});
legend('Location','southeast');
% ylb = 'Median computation time (s)';
% ylabel(ylb);
xticks('')
% yticks([5, 10, 15])

% subplot(212);
nexttile;
hold on;
% grid on;
for i = 2:4
    success_points = median(squeeze(all_error(i, 2, :, :)), 2);
    scatter(all_M, success_points, 48, colors(i, :), 'Marker', marker{i}, 'LineWidth', lineWidth, 'displayname', lgd{i});
    % colors(i, :)
end
% set(gca, 'YScale', 'log')

legend('Location','northeast');
ylb = 'Median computation time (s)';
% label_h = ylabel(ylb, 'Position', [5.7344 0.2300 -1.0000]);
label_h = ylabel(ylb);
label_h.Position(2) = label_h.Position(2) + 0.2;
label_h.Position(1) = label_h.Position(1) - 0.4;
% set(gca, 'ylabelPosition', [100, 100, 100])
xlabel('Number of independent trials')


set(gcf, 'Position',  [100, 100, 700, 300]);
% set(gcf, 'PapserPosition',  [100, 100, 500, 250]);
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
tightfig(gcf)
saveas(gcf, 'Sec_5_2_Four_method_time_2.pdf')
saveas(gcf, '/Users/quanjunlang/Documents/GitHub/Learning-operators-in-Lindblad-Master-Equation-Paper/figures/Sec_5_2_Four_method_time_2.pdf')



%% save data
% save('M_15_79_Nsamples_20.mat')
















%%


% 
% %%
% figure;hold on
% colors = colororder;
% for i = 1:Num_M
%     scatter(all_n, squeeze(log10(all_error(i, :, :))), 500, '.', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :))
% end
% % plot(all_n, all_n.*log(all_n)*trueInfo.rank_RL_true)
% %%
% figure;hold on
% colors = colororder;
% for i = 1:Num_n
%     for k = 1:Num_samples
%         if k==1
%             vs = 'on';
%         else
%             vs = 'off';
%         end
%         lgd = ['n = ', num2str(all_n(i))];
%         scatter(all_M, squeeze(log10(all_error(:, i, k))), 100, 'o', 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), 'DisplayName', lgd,'HandleVisibility', vs)
%     end
% end
% 
% 
% % Convert vector to string for legend
% legendStrings = arrayfun(@num2str, all_n, 'UniformOutput', false);
% % Add the legend using the vector
% % legend(legendStrings, 'Location', 'best');
% legend()
% 
% xlabel('M')
% ylabel('log 10 relative error')
% 
% 
% %%
% figure;hold on;grid on;view(20, 20);
% [XX, YY] = meshgrid(all_n, all_M);
% 
% for k = 1:Num_samples
%     scatter3(XX, YY, squeeze(log10(all_error(:, :, k))), 500, 'k.');
% end
% % scatter3(XX, YY, all_error)
% 
% 
% % Define the grid for x and z
% [x, z] = meshgrid(-2*pi:0.1:2*pi, -2*pi:0.1:2*pi);
% 
% % Calculate y as a function of x
% f = @(x) x.*log(x)*trueInfo.rank_RL_true;
% 
% 
% y = f(XX);
% 
% ee = linspace(log10(min(all_error, [], 'all')), log10(max(all_error, [], 'all')), Num_M);
% [~, ZZ] = meshgrid(all_n, ee);
% surf(XX, y, ZZ)
% ylabel('M')
% xlabel('n')
% 
% 
% 
% %%
% figure;hold on;grid on;view(0, 90);
% all_percent = sum(all_error<0.01, 3)/Num_samples;
% % plot(all_n, f(all_n), 'k', 'LineWidth',20)
% % surf(XX, YY, all_percent)
% 
% % heatmap(all_percent, 'Colormap', jet)
% 
% imagesc(all_n, all_M, all_percent)
% colormap('autumn')
% plot(all_n, f(all_n), 'k', 'LineWidth', 5)
% text(all_n(end-1), f(all_n(end-1))+20, 'M=rnlog(n)')
% 
% % s = 0.2;
% % g = @(x) x.*log(x)*trueInfo.rank_RL_true./x.^s+30;
% % plot(all_n, g(all_n), 'k', 'LineWidth', 5)
% colorbar
% 
% xlabel('System size n')
% ylabel('Trajectory number M')
% title('Recovery rate')
% 
% set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14);
% 
% %% Time analysis
% figure;hold on;grid on;view(20, 20);
% [XX, YY] = meshgrid(all_n, all_M);
% 
% for k = 1:Num_samples
%     scatter3(XX, YY, squeeze(log10(all_time(:, :, k))), 500, 'k.');
% end
% 
% ylabel('M')
% xlabel('n')
% 
% 
% 
