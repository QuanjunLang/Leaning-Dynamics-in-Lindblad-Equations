function [sysInfo, all_rho, trueInfo] = generate_data(sysInfo, varargin)
% Use Python package to generate trajectories


%% Input parser
p = inputParser;
addRequired(p, 'sysInfo');
addOptional(p, 'plotON', 0);

parse(p, sysInfo, varargin{:});
plotON                  = p.Results.plotON;


%%
n = sysInfo.n;
p = sysInfo.p;
tgrid = sysInfo.tgrid;
M = sysInfo.M;
L = sysInfo.L;
if isfield('obs_std', sysInfo)
    obs_std = sysInfo.obs_std;
else
    obs_std = 0;
end

%%
a = pyrunfile("Lindblad.py", 'a', n = n, p = p, tgrid = tgrid, M = M);


%% Load data from Python code
result_temp     = cell(a);
all_traj_temp   = result_temp{1};
H_true          = double(result_temp{2}.full());
C_true_temp     = cell(result_temp{3});
C_true = cell(p, 1);

L_true = -1i*(kron(eye(n), H_true) - kron(H_true.', eye(n)));
E_true = -1i*(vec(eye(n))*vec(H_true.').' - vec(H_true)*vec(eye(n).').');

for i = 1:p
    C_true{i} = double(C_true_temp{i}.full());
    L_C = 0.5*(2*kron(conj(C_true{i}), C_true{i}) - kron((C_true{i}'*C_true{i}).', eye(n)) - kron(eye(n), C_true{i}'*C_true{i}));
    L_true = L_true + L_C;

    E_C = 0.5*(2*vec(C_true{i}')*vec(C_true{i}.').' - vec(C_true{i}'*C_true{i})*vec(eye(n).').' - vec(eye(n).')*vec((C_true{i}'*C_true{i}).').');
    E_true = E_true + E_C;
end





trueInfo.E_true = E_true;
trueInfo.r_true = rank(trueInfo.E_true);
trueInfo.L_true = L_true;
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

trueInfo.all_rho_true = all_rho;
all_rho = all_rho + (randn(size(all_rho)) + 1i*randn(size(all_rho)))*obs_std*abs(mean(all_rho, 'all'));


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

%% plot eigenvalues of E and L
if plotON
    figure;

    subplot(131);
    eigE = eig(trueInfo.E_true);
    plot(eigE, '.', 'markersize', 20);
    xE = max(max(abs(real(eigE))), 1);
    yE = max(max(abs(imag(eigE))), 1);
    xlim([-xE, xE]); ylim([-yE, yE])
    xline(0)
    grid on;
    title('Eigenvalues of E')

    subplot(132);
    eigL = eig(trueInfo.L_true);
    xL = max(max(abs(real(eigL))), 1);
    yL = max(max(abs(imag(eigL))), 1);
    plot(eigL, '.', 'markersize', 20);
    xlim([-xL, xL]); ylim([-yL, yL])
    xline(0)
    grid on;
    title('Eigenvalues of L')


    subplot(133);hold on;
    plot(log10(svd(E_true)), '.-', 'MarkerSize',10);
    plot(log10(svd(L_true)), '.-', 'MarkerSize',10);
    grid on;
end
end


