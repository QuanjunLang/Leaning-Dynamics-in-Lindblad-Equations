function [all_rho, trueInfo, observableInfo] = generate_data(sysInfo, varargin)
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
N_o         = sysInfo.N_o;
tgrid       = sysInfo.tgrid;
FULL_STATE  = sysInfo.FULL_STATE;

%% Generate data using python packages
if FULL_STATE
    a = pyrunfile("Lindblad.py", 'a', n = n, p = p, tgrid = tgrid, M = M);
else
    a = pyrunfile("Lindblad_observable.py", 'a', n = n, p = p, tgrid = tgrid, M = M, N_o = N_o);
end
%% Load data from Python code
result_temp     = cell(a);

% Hamiltonian
H_true          = double(result_temp{2});

% Jump operators
C_true_temp     = cell(result_temp{3});
C_true          = cell(p, 1);
for i = 1:p
    C_true{i} = double(C_true_temp{i});
end

if FULL_STATE
    % Trajectory
    all_rho_temp        = double(result_temp{1});
    all_rho             = permute(all_rho_temp, [3, 4, 2, 1]);
    observableInfo.rho0 = squeeze(all_rho(:, :, 1, :));
else
    % Expectation
    all_expect_temp   = double(result_temp{1});
    all_expect        = permute(all_expect_temp, [2, 3, 1]);
    all_rho           = all_expect;

    % Observable
    O_true_temp     = cell(result_temp{4});
    O_true          = zeros(n, n, N_o, M);
    for m = 1:M
        O_true_temp_m = cell(O_true_temp{m});
        for i = 1:N_o
            O_true(:, :, i, m) = double(O_true_temp_m{i});
        end
    end

    % Initial
    all_initial_temp     = cell(result_temp{5});
    all_initial          = zeros(n, n, M);
    for i = 1:M
        all_initial(:, :, i) = double(all_initial_temp{i});
    end
    
    observableInfo.O    = O_true;
    observableInfo.rho0 = all_initial;
end

trueInfo.computation_time = double(result_temp{end});
%% Construct the matrix of L and the reshaped one E

H_part_true = get_H_part(H_true);

L_true = H_part_true;
E_true = -1i*(vec(eye(n))*vec(H_true.').' - vec(H_true)*vec(eye(n).').');
phi_true = zeros(size(L_true));

for i = 1:p
    L_C = 0.5*(2*kron(conj(C_true{i}), C_true{i}) - kron((C_true{i}'*C_true{i}).', eye(n)) - kron(eye(n), C_true{i}'*C_true{i}));
    L_true = L_true + L_C;

    E_C = 0.5*(2*vec(C_true{i}')*vec(C_true{i}.').' - vec(C_true{i}'*C_true{i})*vec(eye(n).').' - vec(eye(n).')*vec((C_true{i}'*C_true{i}).').');
    E_true = E_true + E_C;

    phi_C = kron(conj(C_true{i}), C_true{i});
    phi_true = phi_true + phi_C;
end



trueInfo.E_true = E_true;
trueInfo.r_true = rank(E_true);
trueInfo.L_true = L_true;
trueInfo.H_true = H_true;
trueInfo.C_true = C_true;
trueInfo.H_part_true = H_part_true;
trueInfo.jump_part_true = L_true - H_part_true;

%% decomposition demo - (Cpt map and kappa)

I = vec(eye(n, n));
kappa_true = 1i*H_true + 0.5*reshape(phi_true'*I, n, n);
kappa_part_true = kron(eye(n, n), kappa_true) + kron(conj(kappa_true), eye(n));

norm(L_true - phi_true + kappa_part_true, 'fro');

trueInfo.phi_true = phi_true;
trueInfo.kappa_true = kappa_true;
trueInfo.kappa_part_true = kappa_part_true;


%% decomposition - (Hamiltonian and Kossakowski)

F = cell(n^2-1, 1);
for l = 1:n-1
    temp = zeros(n, 1);
    temp(1:l) = 1;
    temp(l+1) = -l;
    % F_l = 1i/sqrt(l*(l+1))*(diag(temp));
    F_l = 1/sqrt(l*(l+1))*(diag(temp));
    F{l} = sparse(F_l);
end
l = l+1;
for k = 1:n
    for j = 1:k-1
        e_jk = zeros(n, n);
        e_jk(j, k) = 1;

        F_l = (e_jk + e_jk')/sqrt(2);F{l} = sparse(F_l);l = l+1;
        F_l = -1i*(e_jk - e_jk')/sqrt(2);F{l} = sparse(F_l);l = l+1;
       
    end
end



v = zeros(n^2-1, p);
for q = 1:length(C_true)
    for i = 1:n^2 - 1
        v(i, q) = trace(C_true{q}'*F{i});
    end
end


G = cell(n^2-1, n^2-1);
for k = 1:n^2-1
    for l = 1:n^2-1
        G_kl = kron(conj(F{k}), F{l}) - 0.5*kron((F{k}'*F{l}).', eye(n)) - 0.5*kron(eye(n), F{k}'*F{l});
        G{k, l} = sparse(G_kl);
    end
end


c = zeros(n^2-1, 1);
for j = 1:n^2-1
    c(j) = trace(H_true*F{j}');
end



trueInfo.F = F;
trueInfo.G = G;         % inner product matrix given basis F
trueInfo.c = c;         % coefficient of H on basis F
trueInfo.Kossakowski = v*v';
trueInfo.v_true = v;

%%

% all_rho = zeros(n, n, L, M);
% for m = 1:M
%     % curr_traj_temp = cell(all_traj_temp{m});
%     for t = 1:L
%         % t
%         all_rho(:, :, t, m) = double(all_traj_temp{m}{t}.full());
%         % all_rho(:, :, t, m) = double(curr_traj_temp{t}.full());
%     end
% end

% trueInfo.all_rho_true = all_rho;
% all_rho = all_rho + (randn(size(all_rho)) + 1i*randn(size(all_rho)))*obs_std*abs(mean(all_rho, 'all'));


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

    subplot(221);
    eigE = eig(trueInfo.E_true);
    plot(eigE, '.', 'markersize', 20);
    xE = max(max(abs(real(eigE))), 1);
    yE = max(max(abs(imag(eigE))), 1);
    xlim([-xE, xE]); ylim([-yE, yE])
    xline(0)
    grid on;
    title('Eigenvalues of E')

    subplot(222);
    eigL = eig(trueInfo.L_true);
    xL = max(max(abs(real(eigL))), 1);
    yL = max(max(abs(imag(eigL))), 1);
    plot(eigL, '.', 'markersize', 20);
    xlim([-xL, xL]); ylim([-yL, yL])
    xline(0)
    grid on;
    title('Eigenvalues of L')


    subplot(223);hold on;
    plot(log10(svd(E_true)), '.-', 'MarkerSize',10);
    plot(log10(svd(L_true)), '.-', 'MarkerSize',10);
    grid on;

    subplot(224);hold on;
    plot(exp(eigL), '.', 'markersize', 20);
end
end


