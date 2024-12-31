clc
close all
clear all
random_seed = 0;
rng(random_seed);

addPaths

try
    terminate(pe)
catch
    pause(0.01)
end

% copyright - Quanjun Lang, 2024
%% system settings
sysInfo.n               = 6;           %
sysInfo.M               = 100;          % number of independent trajectories
sysInfo.dt              = 0.000001;     % true data generation time grid
sysInfo.p               = 3;            % number of jumpoperators
sysInfo.steps           = 1;
sysInfo.channel_dt_rate = 1;
sysInfo.observable_option  = 'Channel_first_row_col';

sysInfo.PAPER_FIG_DIR = 'figure';

pe = pyenv(Version='/opt/anaconda3/bin/python', ExecutionMode = 'OutOfProcess');
terminate(pe)


%% Typical case
% sysInfo.M = 80;
% sysInfo = update_sys(sysInfo);
% 
% 
% % Generate data
% if contains(sysInfo.observable_option, 'Channel')
%     % Channel learning
%     [all_rho, trueInfo, observableInfo] = generate_data_channel(sysInfo, 'plotON', 1);
% else
%     % Lindbladian learning
%     [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);
% 
%     b_mat = vec(squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt);
%     b_all = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;
% 
% end
% 
% % all_rho contains the observation data of rho_0 and first-row-col observables
% % This is the setting for ALS-P, ALS-N, ALS-I
% % But not ideal for total ALS, since the sensing matrices are not
% % "incoherent enought"
% 
% % Want to compare with the following four case
% % 1. FFO (Fixed first-row observables):     sUse the same setting with total ALS
% % 2. FRO (Fixed random observables):        Use same initial states, with fixed set of 2N+1 observables
% % 3. CRO (Changing random observables):     Use same initial states, with random 2N+1 observables for each initial state
% % 4. TRO (Total random state observables):  Use totally random initial states and observables, of M(2N+1)-pairs,
% % to match the total amount with previous settings.
% 
% 
% 
% % Extract parameters
% n = sysInfo.n;
% M = sysInfo.M;
% N_o = sysInfo.N_o;
% 
% K = trueInfo.E_true;
% RK = trueInfo.RE_true;
% rank_RK = trueInfo.rank_RE_true;
% RK_sub_blocks = trueInfo.RE_sub_blocks;
% all_rho0    = observableInfo.rho0;
% 
% all_A = cell(4, 1);
% all_b = cell(4, 1);
% 
% 
% [all_A{1}, all_b{1}] = get_A_b_FFO(sysInfo, observableInfo, all_rho);
% [all_A{2}, all_b{2}] = get_A_b_FRO(sysInfo, observableInfo, RK);
% [all_A{3}, all_b{3}] = get_A_b_CRO(sysInfo, observableInfo, RK);
% [all_A{4}, all_b{4}] = get_A_b_TRO(sysInfo, RK);
% 
% compute_dnorm = 0;
% method_list = {'Fixed first row observables', 'Fixed random observables', 'Changed random observables', 'Total random states and observables'};
% 
% 
% RK_est = cell(4, 1);
% RK_Info = cell(4, 1);
% 
% for i = 1:4
%     fprintf('\nJoint ALS for the entire matrix, method = %s...\n', method_list{i})
% 
%     [RK_est{i}, RK_Info{i}] = ALS(all_A{i}, all_b{i}, rank_RK, 'X_true', RK, ...
%         'Nesterov_beta', 0.5, 'loss_tol', 1e-6, 'rel_err_tol', 1e-7, 'debugON', 0, 'maxIter', 500);
% 
%     fprintf('Finished in %.2f seconds \n', RK_Info_total.time)
% 
%     RK_Info{i} = error_analysis(RK_Info{i}, RK_est{i}, RK, 'displayON', 1, 'compute_dnorm', compute_dnorm);
% end




%%

Num_methods = 4;

all_M       = 40:20:80;
Num_M = length(all_M);

Num_samples = 2;


all_RK_Info = cell(Num_M, Num_samples, Num_methods);



for i = 1:Num_M
    for k = 1:Num_samples

        %% Load parameters
        sysInfo.M = all_M(i);
        sysInfo = update_sys(sysInfo);
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        fprintf('The %d-th M = %d out of %d total choice, %d-th sample out of %d\n', i, sysInfo.M, Num_M, k, Num_samples)
        fprintf('\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        %% Generate data
        % Generate data
        if contains(sysInfo.observable_option, 'Channel')
            % Channel learning
            [all_rho, trueInfo, observableInfo] = generate_data_channel(sysInfo, 'plotON', 1);
        else
            % Lindbladian learning
            [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, 'plotON', 1);

            b_mat = vec(squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt);
            b_all = squeeze(all_rho(:, 2, :) - all_rho(:, 1, :))/sysInfo.dt;

        end



        %% Extract parameters
        n = sysInfo.n;
        M = sysInfo.M;
        N_o = sysInfo.N_o;

        K = trueInfo.E_true;
        RK = trueInfo.RE_true;
        rank_RK = trueInfo.rank_RE_true;
        RK_sub_blocks = trueInfo.RE_sub_blocks;
        all_rho0    = observableInfo.rho0;

        all_A = cell(4, 1);
        all_b = cell(4, 1);


        [all_A{1}, all_b{1}] = get_A_b_FFO(sysInfo, observableInfo, all_rho);
        [all_A{2}, all_b{2}] = get_A_b_FRO(sysInfo, observableInfo, RK);
        [all_A{3}, all_b{3}] = get_A_b_CRO(sysInfo, observableInfo, RK);
        [all_A{4}, all_b{4}] = get_A_b_TRO(sysInfo, RK);

        compute_dnorm = 0;
        method_list = {'Fixed first row observables', 'Fixed random observables', 'Changed random observables', 'Total random states and observables'};


        for j = 1:Num_methods
            fprintf('\nJoint ALS for the entire matrix, method = %s...\n', method_list{j})

            [RK_est, all_RK_Info{i, k, j}] = ALS(all_A{i}, all_b{i}, rank_RK, 'X_true', RK, ...
                'Nesterov_beta', 0.5, 'loss_tol', 1e-6, 'rel_err_tol', 1e-7, 'debugON', 0, 'maxIter', 500);

            fprintf('Finished in %.2f seconds \n', all_RK_Info{i, k, j}.time)

            all_RK_Info{i, k, j} = error_analysis(all_RK_Info{i, k, j}, RK_est, RK, 'displayON', 1, 'compute_dnorm', compute_dnorm);
        end

    end
end




%% Extract data
for i = 1:Num_M
    for k = 1:Num_samples
        for j = 1:Num_methods
            all_error(i, k, j) = all_RK_Info{i, k, j}.rel_error_fro_X;
            all_time(i, k, j) = all_RK_Info{i, k, j}.time;
        end
    end
end
%%
hfig = figure;hold on;
colors(1, :) = color1;
colors(2, :) = color2;
colors(3, :) = color3;
colors(4, :) = color4;

marker{1} = '+';
marker{2} = 'o';
marker{3} = '*';
marker{4} = 's';

lineWidth = 2;

lgd{1} = 'FFO';
lgd{2} = 'FRO';
lgd{3} = 'CRO';
lgd{4} = 'TRO';


for i = 1:Num_methods
    success_rate = sum(all_error(:, :, i)<1e-5, 2);
    scatter(all_M, success_rate, 52, 'Color', colors(i, :), 'Marker', marker{i}, 'LineWidth', lineWidth, 'displayname', lgd{i});
end

legend()

%%

figure;hold on;
for i = 1:Num_methods
    mean_time = mean(all_time(:, :, i), 2);
    scatter(all_M, mean_time, 52, 'Color', colors(i, :), 'Marker', marker{i}, 'LineWidth', lineWidth, 'displayname', lgd{i});
end

legend()