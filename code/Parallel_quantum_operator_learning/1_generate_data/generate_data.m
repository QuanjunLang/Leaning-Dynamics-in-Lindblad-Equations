function [all_data, trueInfo, observableInfo] = generate_data(sysInfo, varargin)
% Use Python package to generate trajectories

%% Input parser
p = inputParser;
addRequired(p, 'sysInfo');
addOptional(p, 'plotON', 0);

parse(p, sysInfo, varargin{:});
plotON                  = p.Results.plotON;

%% Load parameters

n           = sysInfo.n;
p           = sysInfo.p;
M           = sysInfo.M;
tgrid       = sysInfo.tgrid;
N_o         = sysInfo.N_o;
observable_option  = sysInfo.observable_option;

%% Generate data using python packages
switch observable_option
    case 'Full_state'
        a = pyrunfile("Lindblad.py", 'a', n = n, p = p, tgrid = tgrid, M = M);
    case 'Multiple_random_observables'
        a = pyrunfile("Lindblad_multiple_observable.py", 'a', n = n, p = p, tgrid = tgrid, M = M, N_o = N_o);
    case 'Single_random_observable'
        a = pyrunfile("Lindblad_single_observable.py", 'a', n = n, p = p, tgrid = tgrid, M = M);
    case 'First_row_col_diag'
        a = pyrunfile("Lindblad_first_row_col_diag.py", 'a', n = n, p = p, tgrid = tgrid, M = M);
end

%% Load data from Python code
result_temp     = cell(a);


all_rho_temp    = double(result_temp{1});
H_true          = double(result_temp{2});   % Hamiltonian
J_true_temp     = cell(result_temp{3});     % Jump operators
J_true          = cell(p, 1);
for i = 1:p
    J_true{i} = double(J_true_temp{i});
end


all_initial_temp     = cell(result_temp{4});  % Initial
all_initial          = zeros(n, n, M);
for i = 1:M
    all_initial(:, :, i) = double(all_initial_temp{i});
end
observableInfo.rho0 = all_initial;


computation_time = double(result_temp{5});


switch observable_option
    case 'Full_state'
        % Trajectory
        all_rho_temp        = double(result_temp{1});
        all_data            = permute(all_rho_temp, [3, 4, 2, 1]);
        observableInfo.rho0 = squeeze(all_data(:, :, 1, :));
    case 'Multiple_random_observables'
        %Expectation
        all_expect_temp   = double(result_temp{1});
        all_data        = permute(all_expect_temp, [2, 3, 1]);

        % Observable
        O_true_temp     = cell(result_temp{6});
        O_true          = zeros(n, n, N_o, M);
        for m = 1:M
            O_true_temp_m = cell(O_true_temp{m});
            for i = 1:N_o
                O_true(:, :, i, m) = double(O_true_temp_m{i});
            end
        end
        observableInfo.O    = O_true;
    case 'Single_random_observable'
        % Expectation
        all_expect_temp   = double(result_temp{1});
        all_data        = permute(all_expect_temp, [2, 1]);

        % Observable
        O_true_temp     = cell(result_temp{6});
        O_true          = zeros(n, n, M);

        for m = 1:M
            O_true(:, :, m) = double(O_true_temp{m});
        end
        observableInfo.O    = O_true;
    case 'First_row_col_diag'
        % Expectation
        all_expect_temp   = double(result_temp{1});
        all_data        = permute(all_expect_temp, [2, 3, 1]);

        % Observable
        O_true_temp     = cell(result_temp{6});
        O_true          = zeros(n, n, N_o);
        for i = 1:N_o
            O_true(:, :, i) = double(O_true_temp{i});
        end
        observableInfo.O    = O_true;
end





% 
% if FULL_STATE
%     % Trajectory
%     all_rho_temp        = double(result_temp{1});
%     all_data             = permute(all_rho_temp, [3, 4, 2, 1]);
%     observableInfo.rho0 = squeeze(all_data(:, :, 1, :));
% else
% 
% 
%     if ~UseDiffObs
% 
%         
%     else
% 

% 
%     end
% 
% 
%     observableInfo.O    = O_true;
% 
% end


%% Construct the matrix of L and the reshaped one E

H_part_true = get_H_part(H_true);

L_true = H_part_true;

for i = 1:p
    L_C = 0.5*(2*kron(conj(J_true{i}), J_true{i}) - kron((J_true{i}'*J_true{i}).', eye(n)) - kron(eye(n), J_true{i}'*J_true{i}));
    L_true = L_true + L_C;
end

trueInfo.L_true = L_true;
trueInfo.RL_true = rearrangement_R(L_true);
trueInfo.rank_RL_true = rank(trueInfo.RL_true);
trueInfo.H_true = H_true;
trueInfo.J_true = J_true;

%% Construct chennel operator
E_true = expm(L_true * sysInfo.channel_dt_rate*sysInfo.dt);

trueInfo.E_true = E_true;
trueInfo.RE_true = rearrangement_R(E_true);
trueInfo.rank_RE_true = rank(trueInfo.RE_true);

%%
RE = trueInfo.RE_true;
RL = trueInfo.RL_true;

n = sysInfo.n;
N_o = sysInfo.N_o;

RE_sub_blocks = cell(n, n);
for i = 1:n
    for j = 1:n
        RE_sub_blocks{i, j} = RE((i-1)*n+1:i*n, (j-1)*n+1:j*n);
    end
end

RE_obs_blocks = cell(N_o, 1);
for i = 1:n
    RE_obs_blocks{i} = RE_sub_blocks{i, i};
end

for j = 1:n-1
    RE_obs_blocks{i + j} = RE_sub_blocks{1, j+1} + RE_sub_blocks{j+1, 1};
end

for k = 1:n-1
    RE_obs_blocks{i + j + k} = 1i*RE_sub_blocks{1, k+1} -1i* RE_sub_blocks{k+1, 1};
end

%%
RL_sub_blocks = cell(n, n);
for i = 1:n
    for j = 1:n
        RL_sub_blocks{i, j} = RL((i-1)*n+1:i*n, (j-1)*n+1:j*n);
    end
end

RL_obs_blocks = cell(N_o, 1);
for i = 1:n
    RL_obs_blocks{i} = RL_sub_blocks{i, i};
end

for j = 1:n-1
    RL_obs_blocks{i + j} = RL_sub_blocks{1, j+1} + RL_sub_blocks{j+1, 1};
end

for k = 1:n-1
    RL_obs_blocks{i + j + k} = 1i*RL_sub_blocks{1, k+1} -1i* RL_sub_blocks{k+1, 1};
end

%%
trueInfo.RE_sub_blocks = RE_sub_blocks;
trueInfo.RL_sub_blocks = RL_sub_blocks;
trueInfo.RE_obs_blocks = RE_obs_blocks;
trueInfo.RL_obs_blocks = RL_obs_blocks;


end





