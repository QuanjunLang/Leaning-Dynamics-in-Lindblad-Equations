clc
clear all
close all
addPaths
% rng(1)

% Multi Matrix ALS test
% Suppose rho_{t+1} = \sum_{k = 1}^r A_k*rho_t*B_k;
% Denote as rho_{t+1} = E(rho_t)
% Estimate {A_k} and {B_k} from multiple observations {rho^{m}_t}
% m = [1, ..., M], t = [1, ..., T]
%
%

%% True para
n = 4;          % Dimension of matrix
r = 3;          % rank of summation
T = 1;         % Total trajectory length
M = 2*r*n;        % Total number of trajectories
tgrid = 1:T+1;  % Time grid

% True parameters A and B, of shape n*n*r
% Where B stores the conjugate transpose of A
A_true = (rand(n, n, r) + 1i*rand(n, n, r))/20;
B_true = zeros(n, n, r);

for i = 1:r
    B_true(:, :, i) = A_true(:, :, i)';
end

% Use K to store the matrix E: n^2 * n^2, maps vec(rho_t) to vec(rho_{t+1})
% E = sum_{k = 1}^r B_k^T \otimes A_k
K_true = zeros(n^2, n^2);
for k = 1:r
    K_true = K_true + kron(B_true(:, :, k).', A_true(:, :, k));
end
% Note that K_true is our goal of estimation, and the unique way of
% decompose K_true as summation of \sum_{k = }^r (V_r')T \otimes V_r
% requires the orthogonality and uniqueness upto a unitary transformation


% Use E to store the matrix \sum_{k = 1}^rvec(A_k.T)vec(B_k), which
% preserve the low rank information
E_true = zeros(n^2, n^2);
for k = 1:r
    E_true = E_true + vec(B_true(:, :, k))*vec(permute(A_true(:, :, k), [2, 1])).';
end
%% Observation data
% Generate Observation data
% Sampling initial condition with Gaussian
all_rho = zeros(n, n, T+1, M);
all_rho(:, :, 1, :) = randn(n, n, 1, M);
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
    figure;hold on;grid on;view(10, 10)
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
% This format helps to analyze the structure of learning A and B from data
% This is essentially a low rank matrix sensing problem, with
% the relationship <P, vec(A.T)vec(B)'>_F = b
% where P is a tensor of n*n*M*T many copies
% The (i, j, t, m)-th matrix of P is
% I_r \otimes 1_{ij} \otimes rho(:,:,t,m)
% each copy contains a matrix of shape rn^2*rn^2
% b is the vector contains n*n*M*T copies of outputs
% and vec(A) and vec(B) are vectors of length rn^2

% test this view point

ind_M = 1;
ind_T = 1;
ind_1 = 2;
ind_2 = 2;

rho0 = all_rho(:, :, ind_T, ind_M);
rho1 = all_rho(ind_1, ind_2, ind_T+1, ind_M);

temp = zeros(n, n);
temp(ind_1, ind_2) = 1;
P = kron(kron(eye(r), temp), rho0);
P_temp = kron(temp, rho0);

% This is the output using rank 1 matrix point of view
vec(permute(A_true, [2, 1, 3])).'*P*vec(B_true)
% This is the correct output
rho1

%% Another way to write the above operations to save time
% A superkron function, help to construct block kroneker product
% Faster and more storage efficient to generate the matrix for regression

PP = superkron(eye(n, n), pagemtimes(A_true, rho0));
PPA = reshape(PP, [n^2. n^2*r]);
PPA*vec(B_true);
vec(rho1);

%% Sanity check: estimate B from A_true (and vice versa)

% tic
% Given A, estimate B
% Construct the regression matrix
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


% right hand side vector b
bB = vec(all_rho(:, :, 2:end, :));

% linear regression, result is nnr*1
% vec_B_est = regmat_B_given_A \ bB;
% reshaping to n*n*r
% toc


% B_est = reshape(vec_B_est, [n, n, r])
% B_true

% Samilar as above, learn A from B_true
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
% A_est = permute(reshape(vec_A_est, [n, n, r]), [2,1,3])
% A_true

% A_est - A_true
% B_est - B_true

%% Things works with mis specified r
true_para.r_true = r;
true_para.A_true = A_true;
true_para.B_true = B_true;
true_para.K_true = K_true;
true_para.E_true = E_true;


r_guess = r+1;
% r
% r_guess
tic
result = multi_mat_als(all_rho, r_guess, true_para);
toc

%%

tic
result = multi_mat_als(all_rho, r_guess);
toc
%% unique decomposition of K into V
E_est = result.E_est;
[U, S, V] = svd(E_est);

%% these are the estimated Kuras matrices
% V_est = result.V_est;
% V1 = V_est(:, :, 1);
% V2 = V_est(:, :, 2);
% V3 = V_est(:, :, 3);
%
% trace(V1*V2')
% trace(V1*V1')
% Which are orthogonal

%% Singular value thresholding
% Seems to be a faster algorithm. Let's see






ind_M = 1;
ind_T = 1;
ind_1 = 2;
ind_2 = 2;

rho0 = all_rho(:, :, ind_T, ind_M);
rho1 = all_rho(ind_1, ind_2, ind_T+1, ind_M);

temp = zeros(n, n);
temp(ind_1, ind_2) = 1;
P = kron(kron(eye(r), temp), rho0);

% This is the output using rank 1 matrix point of view
vec(permute(A_true, [2, 1, 3])).'*P*vec(B_true)          
% This is the correct output
rho1

all_P = zeros(n^2, n^2, n, n, T, M);
for t = 1:T
    for m = 1:M
        rho0 = all_rho(:, :, t, m);
        for i = 1:n
            for j = 1:n
                temp = zeros(n, n);
                temp(i, j) = 1;
                all_P(:, :, i, j, t, m) = kron(temp, rho0);
            end 
        end
    end
end

%%
A = reshape(all_P, [n^4, n*n*M*T])';
X = reshape(E_true, [], 1);
b = conj(reshape(squeeze(all_rho(:, :, 2:end, :)), [], 1));



%%
tic

tau = 0.05;
delta = 0.0005;
niter = 100;

rk = zeros(niter, 1);
err = zeros(niter, 1);

all_X = zeros(n^2, n^2, niter);
all_Y = zeros(n^2*M*T, 1, niter);

all_X(:, :, 1) = rand(n^2, n^2);
% all_Y(:, :, 1) = rand(n^2*M*T, 1);

At = A';
for k = 1:niter
    X0 = all_X(:, :, k);
    Y0 = all_Y(:, :, k);
    
    X1_vec  = At*Y0;
    X1_temp = reshape(X1_vec, n^2, n^2);
    % [U, S, V] = svt(X1_temp, 'lambda' ,tau);

    [U, S, V] = svd(X1_temp);
    % diag(S)
    S(S<tau) = 0;
    % diag(S)
    X1 = U*S*V';

    X1_vec = reshape(X1, [], 1);
    Y1 = Y0 + delta*(b - A*X1_vec);
    
    all_X(:, :, k+1) = X1;
    all_Y(:, :, k+1) = Y1;

    err(k) = norm(A*X1_vec - b);
    rk(k) = sum(S>0, 'all');
end

toc

figure;hold on;
plot(log10(err))
plot(rk)




