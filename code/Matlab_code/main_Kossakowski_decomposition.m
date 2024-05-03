clc
close all
clear all
rng(1)
addPaths

%%
sysInfo.n       = 3;            %
sysInfo.M       = 1;           % number of independent trajectories
sysInfo.dt      = 0.05;         % true data generation time grid
sysInfo.p       = 2;            % number of jump operators


sysInfo.steps   = 500;
sysInfo = update_sys(sysInfo);

[all_rho, trueInfo] = generate_data(sysInfo);


obsInfo.obs_std = 1e-4;
obsInfo.obs_gap = 20;
obsInfo.obs_len = 10;


[all_rho_obs, obsInfo] = generate_observation_data(all_rho, sysInfo, obsInfo);


%%
L = trueInfo.L_true;

Kossakowski = trueInfo.Kossakowski;
Hamiltonian = trueInfo.H_true;
H_part = trueInfo.H_part_true;
v_true = trueInfo.v_true;
p = sysInfo.p;
F = trueInfo.F;
n = sysInfo.n;

G = cell(n^2-1, n^2-1);

for k = 1:n^2-1
    for l = 1:n^2-1
        G_kl = kron(conj(F{k}), F{l}) - 0.5*kron((F{k}'*F{l}).', eye(n)) - 0.5*kron(eye(n), F{k}'*F{l});
        G{k, l} = sparse(G_kl);
    end
end

%% iterations
clc
result = L_decomposition_hamiltonian_kossakowski(L, p, trueInfo);











%%


K_part = zeros(size(H_part));

for k = 1:n^2 - 1
    for l = 1:n^2 - 1
        K_part = K_part + Kossakowski(k, l)*G{k, l};
    end
end

norm(K_part + H_part - L);

%% test thte vectorized equivalence
G_v = cell(n^2-1, p);
v = v_true;

for k = 1:n^2-1
    for q = 1:p
        temp = sparse(zeros(n^4, 1));
        for l = 1:n^2 - 1
            temp = temp + conj(v(l, q))*vec(G{k, l});
        end
        G_v{k, q} = temp;
    end
end



K_part_vv = zeros(n^4, 1);

for k = 1:n^2-1
    for q = 1:p
        K_part_vv = K_part_vv + v(k, q)*G_v{k, q};
    end
end

H_part_vec = vec(H_part);
L_vec = vec(L);


norm(L_vec - H_part_vec - K_part_vv)



%% vectorize u

% 
% 
% A_u = zeros(n^4, (n^2-1)*p);
% u = v_true;
% v = v_true';
% 
% s = 1;
% for q = 1:p
%     for k = 1:n^2-1
% 
%         temp = sparse(zeros(n^4, 1));
%         for l = 1:n^2 - 1
%             temp = temp + v(q, l)*vec(G{k, l});
%         end
%         A_u(:, s) = temp;s = s+1;
%     end
% end
% 
% 
% norm(A_u\(L_vec - H_part_vec) - vec(u))
v = v_true;
u = v_true;
u_est = L_decomposition_Kossakowski_u(L_vec, v, Hamiltonian, G, n, p);
norm(u_est - u)
%% vectorize v
A_v = zeros(n^4, (n^2-1)*p);
u = v_true;
v = v_true;

s = 1;
for q = 1:p
    for l = 1:n^2-1
        temp = sparse(zeros(n^4, 1));
        for k = 1:n^2 - 1
            temp = temp + u(k, q)*vec(G{k, l});
        end
        A_v(:, s) = temp;s = s+1;
    end
end


% norm(conj(A_v\(L_vec - H_part_vec)) - vec(v))
v_est = L_decomposition_Kossakowski_v(L_vec, u, Hamiltonian, G, n, p);
norm(v_est - v)
%%




E = zeros(n^2- 1, n^2-1);
for i = 1:n^2-1
    for j = 1:n^2-1
        E(i, j) = trace(F{i}*F{j}');
    end
end


%% Learning Hamiltonian
C = trueInfo.C_true;
H = trueInfo.H_true;
L_true = trueInfo.L_true;
% all_rho = zeros(n, n, 2, n^2);
% l = 1;
% for i = 1:n
%     for j = 1:n
%         temp = sparse(zeros(n, n));
%         temp(i, j) = 1;
%         all_rho(:, :, 1, l) = temp;
%         all_rho(:, :, 1, l) = randn(n, n);
%         all_rho(:, :, 2, l) = reshape(L_true*vec(temp), [n, n]);
%         l = l+1;
%     end
% end
% M = n^2;

M = 1;
for m = 1:M
    % all_rho(:, :, 1, m) = temp;
    all_rho(:, :, 1, m) = randn(n, n);
    all_rho(:, :, 2, m) = reshape(L_true*vec(all_rho(:, :, 1, m)), [n, n]);
    % l = l+1;
end
% end



% M = n^2;
L = 2;
use_Lindbladian = 1;


all_A = cell(L-1, M);
all_b = cell(L-1, M);

for m = 1:M
    for l = 1:L-1
        curr_rho = all_rho(:, :, l, m);
        next_rho = all_rho(:, :, l+1, m);
        all_A{l, m} = kron(curr_rho.', eye(n, n)) - kron(eye(n, n), curr_rho);

        if use_Lindbladian
            Lindbladian = combine_Lindbladian(C, curr_rho);
            b_temp = 1i*(next_rho - Lindbladian);
        else
            b_temp = next_rho * 1i;
        end
        all_b{l, m} = reshape(b_temp, [], 1);
    end
end

AA = cat(1, all_A{:});
bb = cat(1, all_b{:});


H_est = reshape(AA\bb, [n, n]);

norm(AA*vec(H) - bb)
norm(AA*vec(H_est) - bb)
H - H_est
%%
% clc
L_true = trueInfo.L_true;
rho = curr_rho;

mat = kron(curr_rho.', eye(n, n)) - kron(eye(n, n), curr_rho);

H_part = -1i*mat*vec(H);
L_part = L_true*vec(curr_rho);
C_part = vec(Lindbladian);

norm(L_part - H_part - C_part)



H_rho = -1i*(H*rho - rho*H);
L_rho = reshape(L_true*vec(curr_rho), [n, n]);
C_rho = Lindbladian;


norm(L_rho - H_rho - C_rho)

%% use double to store G. Too expensive


%
% tic
% G = zeros(n^2, n^2, n^2-1, n^2-1);
%
% for i = 1:n^2-1
%     for j = 1:n^2-1
%         G(:, :, i, j) = kron(conj(F{i}), F{j}) - 0.5*kron((F{i}'*F{j}).', eye(n)) - 0.5*kron(eye(n), F{i}'*F{j});
%     end
% end
%
%
% K_part = zeros(size(H_part));
%
% for k = 1:n^2 - 1
%     for l = 1:n^2 - 1
%         K_part = K_part + Kossakowski(k, l)*G(:, :, k, l);
%     end
% end
%
% norm(K_part + H_part - L)
%

%% Learning Hamiltonian
C = trueInfo.C_true;
H = trueInfo.H_true;
L_true = trueInfo.L_true;
% all_rho = zeros(n, n, 2, n^2);
% l = 1;
% for i = 1:n
%     for j = 1:n
%         temp = sparse(zeros(n, n));
%         temp(i, j) = 1;
%         all_rho(:, :, 1, l) = temp;
%         all_rho(:, :, 1, l) = randn(n, n);
%         all_rho(:, :, 2, l) = reshape(L_true*vec(temp), [n, n]);
%         l = l+1;
%     end
% end
% M = n^2;

M = 1;
for m = 1:M
    % all_rho(:, :, 1, m) = temp;
    all_rho(:, :, 1, m) = randn(n, n);
    all_rho(:, :, 2, m) = reshape(L_true*vec(all_rho(:, :, 1, m)), [n, n]);
    % l = l+1;
end
% end



% M = n^2;
L = 2;
use_Lindbladian = 1;


all_A = cell(L-1, M);
all_b = cell(L-1, M);

all_real_A = cell(L-1, M);
all_real_b = cell(L-1, M);

c = zeros(n^2-1, 1);
for j = 1:n^2-1
    c(j) = trace(H*F{j}');
end

H_proj = zeros(size(H));
for j = 1:n^2-1
    H_proj = H_proj + c(j)*F{j};
end


for m = 1:M
    for l = 1:L-1
        curr_rho = all_rho(:, :, l, m);
        next_rho = all_rho(:, :, l+1, m);
        all_A{l, m} = kron(curr_rho.', eye(n, n)) - kron(eye(n, n), curr_rho);
        
        temp = zeros(n^2, n^2-1);
        for j = 1:n^2 - 1
            temp(:, j) = all_A{l, m}*vec(F{j});
        end
        A_R = real(temp);
        A_I = imag(temp);
        A = [A_R; A_I];
        all_real_A{l, m} = A;


        if use_Lindbladian
            Lindbladian = combine_Lindbladian(C, curr_rho);
            b_temp = 1i*(next_rho - Lindbladian);
        else
            b_temp = next_rho * 1i;
        end
        all_b{l, m} = reshape(b_temp, [], 1);

        b_R = real(all_b{l, m});
        b_I = imag(all_b{l, m});
        b = [b_R; b_I];
        all_real_b{l, m} = b;
    end
end

AA = cat(1, all_A{:});
bb = cat(1, all_b{:});

real_AA = cat(1, all_real_A{:});
real_bb = cat(1, all_real_b{:});


% H_est = reshape(AA\bb, [n, n]);

% norm(AA*vec(H) - bb)
% norm(AA*vec(H_est) - bb)
% H - H_est


c_est = real_AA \ real_bb;
norm(c - c_est)

