clc
close all
clear all
rng(0)
addPaths
pe = pyenv(Version='/usr/bin/python3', ExecutionMode = 'OutOfProcess');
terminate(pe)

%% Rank of E is 2+N_J


%%


all_n = 2:1:10;
all_p = 1:1:10;
all_r = zeros(length(all_n), length(all_p));
for i = 1:length(all_n)
    for j = 1:length(all_p)
        fprintf('i=%d, j=%d\n', i, j)
        sysInfo.n       = all_n(i);            %
        sysInfo.M       = 1;            % number of independent trajectories
        sysInfo.dt      = 0.1;        % true data generation time grid
        sysInfo.p       = all_p(j);              % number of jump operators
        sysInfo.steps   = 10;
        sysInfo = update_sys(sysInfo);
        [all_rho, trueInfo] = generate_data(sysInfo);

        all_r(i, j) = trueInfo.r_true;
    end
end







