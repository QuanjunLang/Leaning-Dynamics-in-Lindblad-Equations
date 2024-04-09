clc
close all
clear all
rng(1)
% pyenv(Version='/opt/anaconda3/envs/py3.11/bin/python')


sysInfo.n       = 5;
sysInfo.steps   = 1;
sysInfo.dt      = 0.0001;
sysInfo.p       = 0;
sysInfo.M       = 1000;

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

%% Demo Learning Hermitian

% Generate data
[sysInfo, all_rho, trueInfo] = generate_data(sysInfo);

% Learning Hermitian
H_est = learning_Hermitian(sysInfo, all_rho);
H_err = Hermitian_error(H_est, trueInfo.H_true)

H_est = learning_Hermitian(sysInfo, all_rho, trueInfo.C_true);
H_err = Hermitian_error(H_est, trueInfo.H_true)

a = 1;
% %% Demo Learnng Lindbladian
% 
% 
% ind = 1;
% C_true = trueInfo.C_true;
% H_true = trueInfo.H_true;
% 
% 
% C_initial = C_true{1};
% C_initial = rand(sysInfo.n, sysInfo.n);
% 
% C_est = learning_Lindbladian_ind(sysInfo, all_rho, trueInfo.H_true, C_initial, {}, C_true);
% 
% 
% 









%% Loop to check performance
N_dt    = 10;
N_n     = 5;

all_dt  = 10.^linspace(-7, -1, N_dt);
all_n   = floor(10.^linspace(0.4, 1.7, N_n));


all_err = zeros(N_dt, N_n);
all_time = zeros(N_dt, N_n);
%%
for i = 1:N_dt
    for j = 1:N_n
        [i ,j]
        sysInfo.n = all_n(j);
        sysInfo.dt = all_dt(i);
        sysInfo = update_sys(sysInfo);

        [sysInfo, all_rho, trueInfo] = generate_data(sysInfo);

        tic;
        H_est = learning_Hermitian(sysInfo, all_rho);
        H_err = Hermitian_error(H_est, trueInfo.H_true);

        all_time(i, j) = toc;
        all_err(i, j) = H_err;
    end
end


%%
plot(log10(all_dt), log10(all_err))
legend()



%%
figure
plot(log10(all_n), log10(all_time))
legend()

%% Learning Lindbladian












% %% sanity check
% % A       = all_A{1, 1};
% rho     = all_rho(:, :, 1, 1);
% d_rho   = (all_rho(:, :, 2, 1) - all_rho(:, :, 1, 1))/dt
% H = H_true;
% 
% increment_mat = -(H*rho - rho*H)*1i
% -(H*rho - rho*H)*1i - d_rho;
% 
% %% sanity check
% A = kron(rho.', eye(n, n)) - kron(eye(n, n), rho)
% vec_H = reshape(H_true, [], 1)
% 
% A * vec_H *(-1i)
% reshape(d_rho, [], 1)
% 
% 
% %%
% A1 = kron(rho.', eye(n, n))
% A1*vec_H_2
% 
% reshape(H_true*rho, [], 1)
%%
% H_true  = rand
% H_true  = [3, 3-2i;3+2i, 2];

% Lindbladian_ON = 0;
% if Lindbladian_ON
    
    % C_true  = zeros(n, n, p);
    % C_true(:, :, 1) = 0.1*[1, 0;0, -2];
    % C_true(:, :, 2) = 0.2*[0, -i;i, 0];
% else
%     C_true = [];
% end

% rho0    = [0.25, 2+i;2-i, 0.75];
% rho0    = rherm(n);
% rho0    = rho0/trace(rho0);
% a = pyrunfile("Lindblad.py", 'a', h = H_true, tgrid = tgrid, rho0 = rho0, c = C_true);











