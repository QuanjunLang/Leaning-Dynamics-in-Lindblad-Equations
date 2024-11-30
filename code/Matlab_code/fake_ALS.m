r = trueInfo.r_true;
E = trueInfo.E_true;

[U_temp, S, V_temp] = svd(E);
UU = U_temp(:, 1:r)*sqrt(S(1:r, 1:r));
VV = V_temp(:, 1:r)*sqrt(S(1:r, 1:r));



%%
all_rho = all_pair_data.prony;

[n, ~, T_temp, M] = size(all_rho);
T = T_temp - 1;

% % Given Data
% m = 10; % Number of matrices A_i and scalars b_i
% n = 5;  % Dimension of U
% p = 5;  % Dimension of V
% r = 3;  % Rank of the factorization

% Example input data (random initialization)
A = zeros(n^2, n^2, n^2*M*T);  % A_i matrices
for i = 1:n
    for j = 1:n
        E_ij = zeros(n, n);
        E_ij(i, j) = 1;
        A_temp = superkron(E_ij, squeeze(all_rho(:, :, 1, :)));
        ind = ((i-1)*n+j-1)*M*T+1:((i-1)*n+j)*M*T;
        A(:, :, ind) = A_temp;
    end
end

b = zeros(n^2*M*T, 1);
for i = 1:n
    for j = 1:n
        b_temp = squeeze(all_rho(i, j, 2, :));
        ind = ((i-1)*n+j-1)*M*T+1:((i-1)*n+j)*M*T;
        b(ind) = b_temp;
    end
end

m = length(A);


%%




loss_function(UU, VV, A, b, m)




VV_opt = solve_V_least_squares(UU, A, b, m, n^2, r)
VV


UU_opt = solve_U_least_squares(VV, A, b, m, n^2, r)
UU


%%

% % Generate true U and V
% true_U = rand(n, r);  % True U
% true_V = rand(p, r);  % True V
% 
% % Compute true b without noise
% true_b = zeros(m, 1);
% for i = 1:m
%     true_b(i) = trace(true_U * true_V' * A(:,:,i));
% end
% 
% % Add Gaussian noise to b
% noise_std = 0.1;  % Standard deviation of the noise
% b = true_b + noise_std * randn(m, 1);

% Display the noisy b
disp('Noisy b:');
disp(b);

% Initialization of U and V for optimization
U = rand(n^2, r) + 1i*rand(n^2, r);  % Initial U
V = rand(n^2, r) + 1i*rand(n^2, r);  % Initial V

U = UU;
V = VV;

% % Initialization of complex U and V for optimization
% U = rand(n, r) + 1i * rand(n, r);  % Initial complex U
% V = rand(p, r) + 1i * rand(p, r);  % Initial complex V

% Parameters for the optimization
max_iter = 1000;  % Maximum number of iterations
tol = 1e-10;       % Convergence tolerance
all_loss = zeros(max_iter, 1);

% Alternating minimization loop using fminunc
options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'off', 'SpecifyObjectiveGradient', false);

objective = @(U, V) loss_function(U, V, A, b, m);

for iter = 1:max_iter
    iter
    % Step 1: Fix V, optimize U using fminunc
    % U_vec0 = [real(U(:)); imag(U(:))]; % Concatenate real and imaginary parts
    % U_vec_opt = fminunc(@(U_vec) loss_U(U_vec, V, A, b, m, n, r), U_vec0, options);
    % U = reshape(U_vec_opt(1:end/2), n^2, r) + 1i * reshape(U_vec_opt(end/2+1:end), n^2, r);  % Reconstruct complex U

    U = solve_U_least_squares(V, A, b, m, n^2, r);
    
    % Step 2: Fix U, optimize V using fminunc
    % V_vec0 = [real(V(:)); imag(V(:))]; % Concatenate real and imaginary parts
    % V_vec_opt = fminunc(@(V_vec) loss_V(U, V_vec, A, b, m, n, r), V_vec0, options);
    % V = reshape(V_vec_opt(1:end/2), n^2, r) + 1i * reshape(V_vec_opt(end/2+1:end), n^2, r);  % Reconstruct complex V
    % 
    V = solve_V_least_squares(U, A, b, m, n^2, r);
    % Check for convergence
    if iter > 1 && abs(objective(U, V) - prev_loss) < tol
        break;
    end
    prev_loss = objective(U, V);
    all_loss(iter) = prev_loss;
end

% Display results
disp('Optimized U:');
disp(U);
disp('Optimized V:');
disp(V);
disp('Final loss value:');
disp(objective(U, V));

%% Loss Function Definitions
function loss = loss_function(U, V, A, b, m)
    % Compute the first term (least squares)
    term1 = 0;
    for i = 1:m
        term1 = term1 + abs(trace(U * V' * A(:,:,i)) - b(i))^2;
    end
    term1 = term1 / (2 * m);
    
    % Compute the second term (regularization)
    term2 = norm(U' * U - V' * V, 'fro')^2 / 8;
    
    % Total loss
    loss = term1 + term2;
end

% Loss function for U (with V fixed)
function loss_U_val = loss_U(U_vec, V, A, b, m, n, r)
    % Separate real and imaginary parts
    U_real = reshape(U_vec(1:end/2), n^2, r);
    U_imag = reshape(U_vec(end/2+1:end), n^2, r);
    U = U_real + 1i * U_imag;  % Reconstruct complex U
    
    loss_U_val = loss_function(U, V, A, b, m);  % Use the full loss function with fixed V
end

% Loss function for V (with U fixed)
function loss_V_val = loss_V(U, V_vec, A, b, m, n, r)
    % Separate real and imaginary parts
    V_real = reshape(V_vec(1:end/2), n^2, r);
    V_imag = reshape(V_vec(end/2+1:end), n^2, r);
    V = V_real + 1i * V_imag;  % Reconstruct complex V
    
    loss_V_val = loss_function(U, V, A, b, m);  % Use the full loss function with fixed U
end


%% Solve for U using least squares with V fixed
function U_opt = solve_U_least_squares(V, A, b, m, n, r)
    % Stack the system for least squares solution
    A_stack = zeros(n*r, m);
    b_stack = zeros(m, 1);
    
    for i = 1:m
        A_stack(:, i) = vec(A(:,:,i)'*V);
        b_stack(i) = b(i);
    end
    
    % Solve the least squares problem
    U_vec = A_stack' \ b_stack;  % Least squares solution
    U_opt = reshape(U_vec, n, r); % Reshape back to matrix form
end

%% Solve for V using least squares with U fixed
function V_opt = solve_V_least_squares(U, A, b, m, n, r)
    % Stack the system for least squares solution
    A_stack = zeros(n*r, m);
    b_stack = zeros(m, 1);
    
    for i = 1:m
        A_stack(:, i) = vec(A(:,:,i)' * U);
        b_stack(i) = b(i);
    end
    
    % Solve the least squares problem
    V_vec = A_stack' \ b_stack;  % Least squares solution
    V_opt = reshape(V_vec, n, r); % Reshape back to matrix form
end


%% Vectorize a matrix into a column vector
function vecA = vec(A)
    vecA = A(:);
end