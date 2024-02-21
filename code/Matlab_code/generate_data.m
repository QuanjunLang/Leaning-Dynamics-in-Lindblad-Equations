function [sysInfo, all_rho, trueInfo] = generate_data(sysInfo)
% Use Python package to generate trajectories

%%
n = sysInfo.n;
p = sysInfo.p;
dt = sysInfo.dt;
tgrid = sysInfo.tgrid;
M = sysInfo.M;
L = sysInfo.L;


%%
a = pyrunfile("Lindblad.py", 'a', n = n, p = p, tgrid = tgrid, M = M);


%% Load data from Python code
result_temp     = cell(a);
all_traj_temp   = result_temp{1};
H_true          = double(result_temp{2}.full());
C_true_temp     = cell(result_temp{3});
C_true = cell(p, 1);
for i = 1:p
    C_true{i} = double(C_true_temp{i}.full());
end

trueInfo.H_true = H_true;
trueInfo.C_true = C_true;
%%
all_rho = zeros(n, n, L, M);
for m = 1:M
    curr_traj_temp = cell(all_traj_temp{m});
    for t = 1:L
        all_rho(:, :, t, m) = double(curr_traj_temp{t}.full());
    end
end

%% plot sample trajectories

% if plotON
%     figure;
%     hold on;
%     grid on;
%     view(10, 10)
% 
%     rho = all_rho(:, :, :, 1);
%     for i = 1:n
%         for j = 1:n
%             plot3(real(squeeze(rho(i, j, :))), imag(squeeze(rho(i, j, :))), tgrid, 'linewidth', 3, ...
%                 'DisplayName',['(', num2str(i), ',' num2str(j), ')']);
% 
%         end
%     end
%     legend()
% 
% end




end


