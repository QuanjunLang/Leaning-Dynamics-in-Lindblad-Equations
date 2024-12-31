function [all_data, trueInfo, observableInfo] = generate_data_channel(sysInfo, varargin)
% Use Python package to generate trajectories

total_tic = tic;
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

random_seed = randi(10000);
%% Generate data using python packages
fprintf('Data generation,  with n = %d, M = %d and %d Kraus operators\n', n, M, p);
switch observable_option
    case 'Channel_first_row_col'
        a = pyrunfile("Channel_first_row_col.py", 'a', n = n, p = p, M = M, random_seed=random_seed);
end

%% Load data from Python code
result_temp     = cell(a);
trueInfo.computation_time = double(result_temp{5});
fprintf('Data generation finished in %.2f seconds\n', trueInfo.computation_time);



E_true          = double(result_temp{2});   % Channel operator
Choi_true       = double(result_temp{3});     % Jump operators

all_rho1        = permute(double(result_temp{7}), [2, 3, 1]);


Kraus_true_temp = cell(result_temp{8});     % Jump operators
Kraus_true      = cell(p, 1);
for i = 1:p
    Kraus_true{i} = double(Kraus_true_temp{i});
end


try all_initial_temp     = cell(result_temp{4});  % Initial
    all_initial          = zeros(n, n, M);
    for i = 1:M
        all_initial(:, :, i) = double(all_initial_temp{i});
    end
catch ME
    % Handle the error
    % disp('An error occurred during the conversion:');
    all_initial_temp = double(result_temp{4});
    all_initial = permute(all_initial_temp, [2, 3, 1]);
end

observableInfo.rho0 = all_initial;


switch observable_option
    case {'Channel_first_row_col'}
        % Expectation
        all_expect_temp   = double(result_temp{1});
        all_data        = permute(all_expect_temp, [2, 1]);

        % Observable
        O_true_temp     = cell(result_temp{6});
        O_true          = zeros(n, n, N_o);
        for i = 1:N_o
            O_true(:, :, i) = double(O_true_temp{i});
        end
        observableInfo.O        = O_true;
        observableInfo.all_rho1 = all_rho1;
end



%% Construct the matrix of L and the reshaped one E

trueInfo.E_true = E_true;
trueInfo.RE_true = rearrangement_R(E_true);
trueInfo.rank_RE_true = rank(trueInfo.RE_true);
trueInfo.Kraus_true = Kraus_true;
trueInfo.Choi_true = Choi_true;
%%
RE = trueInfo.RE_true;

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
trueInfo.RE_sub_blocks = RE_sub_blocks;
trueInfo.RE_obs_blocks = RE_obs_blocks;

fprintf('Data generation total time: %.2f seconds \n', toc(total_tic));
end





