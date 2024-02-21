clc
clear all
close all
addPaths
% rng(1)

% Multi Matrix ALS test
% Suppose rho_{t+1} = \sum_{k = 1}^r A_j*rho_t*B_k;
% Estimate {A_k} and {B_k} from observations {rho_t}

%% True para
n = 2;          % Dimension of matrix
r = 2;          % rank of summation 
T = 2;        % Total trajectory length
M = 2;          % Total number of trajectories
tgrid = 1:T+1;  % Time grid

% True para
A_true = (rand(n, n, r) + 1i*rand(n, n, r))/20;
B_true = zeros(n, n, r);

for i = 1:r
    B_true(:, :, i) = A_true(:, :, i)';
end

% Observation data
all_rho = zeros(n, n, T+1, M);
all_rho(:, :, 1, :) = rand(n, n, 1, M);
for t = 1:T
    all_rho0 = all_rho(:, :, t, :);
    all_rho1 = zeros(n, n, 1, M);
    for k = 1:r
        temp = pagemtimes(A_true(:, :, k), all_rho0);
        all_rho1 = all_rho1 + pagemtimes(temp, B_true(:, :, k));
    end
    all_rho(:, :, t+1, :) = all_rho1;
end

%% Plot trajectory

plotON = 1;
if plotON
    figure;
    hold on;
    grid on;
    view(10, 10)

    rho = all_rho(:, :, :, 1);
    for i = 1:n
        for j = 1:n
            plot3(real(squeeze(rho(i, j, :))), imag(squeeze(rho(i, j, :))), tgrid, 'linewidth', 3, ...
                'DisplayName',['(', num2str(i), ',' num2str(j), ')']);

        end
    end
    legend()

end

%% Sanity Check: Vectorize the matrix A and B and write it as matrix multiplication
rho0 = all_rho(:, :, 1, 1);
rho1 = all_rho(:, :, 2, 1)
temp = zeros(n, n);
temp(2, 2) = 1;
P = kron(kron(eye(r), temp), rho0)

% vec(A_true(:, :, 1))'*P * vec(B_true(:, :, 1))  + vec(A_true(:, :, 2))'*P * vec(B_true(:, :, 2))
vec(permute(A_true, [2, 1, 3])).'*P*vec(B_true)
% vec(A_true(:, :, 1))'*P * vec(B_true(:, :, 1))  + vec(A_true(:, :, 2))'*P * vec(B_true(:, :, 2))

%% Another way to write the above operations to save time

PP = superkron(eye(n, n), pagemtimes(A_true, rho0));
PPA = reshape(PP, [n^2. n^2*r]);
PPA*vec(B_true)

vec(rho1)




%% Compare two ways to construct the regression matrix

tic
regmat_B_given_A = zeros(M*T*n^2, r*n^2);
k = 1;
for m = 1:M
    for t = 1:T
        rho0 = all_rho(:, :, t, m);
        PP = superkron(eye(n, n), pagemtimes(A_true, rho0));
        PPA = reshape(PP, [n^2. n^2*r]);
        regmat_B_given_A((k-1)*n^2+1:k*n^2, :) = PPA;
        k = k+1;
    end
end
toc

bB = vec(all_rho(:, :, 2:end, :));

vec_B_est = regmat_B_given_A \ bB;
B_est = reshape(vec_B_est, [n, n, r])
B_true


tic
regmat_A_given_B = zeros(M*T*n^2, r*n^2);
k = 1;
for m = 1:M
    for t = 1:T
        rho0 = all_rho(:, :, t, m);
        PP = superkron(eye(n, n), pagemtimes(permute(B_true, [2, 1, 3]), permute(rho0, [2,1,3,4])));
        PPA = reshape(PP, [n^2. n^2*r]);
        regmat_A_given_B((k-1)*n^2+1:k*n^2, :) = PPA;
        k = k+1;
    end
end
toc

bA = vec(permute(all_rho(:, :, 2:end, :), [2,1,3,4]));

vec_A_est = regmat_A_given_B \ bA;
A_est = permute(reshape(vec_A_est, [n, n, r]), [2,1,3])
A_true

A_est - A_true
B_est - B_true
%% Iteration

% build tensor X
% X = zeros(0

niter = 20;

all_A = zeros(r*n^2, niter);    % u = (vec(A_1^T)^T, ..., vec(A_r^T)^T))
all_B = zeros(r*n^2, niter);    % v = (vec(B_1)^T, ..., vec(B_r)^T))



for i = 1:niter
    % Estimate A first 




    AA_U = cell(T, 1);
    bb = cell(T, 1);

    for k = 1:T
        AA_U{k} = kron((AA(:, :, k)*V_cur).', eye(m, m));
        bb{k} = reshape(all_rho(:, :, k), [], 1);
    end
    A_U = cat(1, AA_U{:});
    b_U = cat(1, bb{:});

    U_est = reshape(A_U\b_U, [m, m]);
    % U_est = U_est./U_est(1, 1);

    all_vec_A(:, :, i) = U_est;
    U_cur = U_est;

    % % Estimate V next 
    AA_V = cell(T, 1);
    % bb = cell(K, 1);

    for k = 1:T
        AA_V{k} = kron(eye(n, n), U_cur*AA(:, :, k));
        % bb{k} = reshape(B(:, :, k), [], 1);
    end

    A_V = cat(1, AA_V{:});
    b_V = cat(1, bb{:});

    V_est = reshape(A_V\b_V, [n, n]); 

    % V_est = U_est;
    %         \\\\\\\\\\\\
    all_vec_B(:, :, i) = V_est;
    V_cur = V_est;


    err_B(i) = norm(all_rho - pagemtimes(pagemtimes(U_est, AA), V_est), 'fro');
end


err_V = squeeze(sum((all_vec_B - B_true).^2, [1,2]));
err_U = squeeze(sum((all_vec_A - A_true).^2, [1,2]));


figure;hold on
plot(log10(err_V), ':', 'LineWidth', 3)
plot(log10(err_U))
plot(log10(err_B))

%%
U_est = all_vec_A(:,  :, end);
V_est = all_vec_B(:,  :, end);


B_est = pagemtimes(pagemtimes(U_est, AA), V_est);

% err_B = 

%%


