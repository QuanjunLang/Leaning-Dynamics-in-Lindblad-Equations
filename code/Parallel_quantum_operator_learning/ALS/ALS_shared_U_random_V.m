function [M_est, outputInfo] = ALS_shared_U_random_V(A, all_b, r, varargin)

%%
p = inputParser;

addRequired(p, 'A');
addRequired(p, 'all_b');
addRequired(p, 'r');
addOptional(p, 'debugON', 0); 
addOptional(p, 'plotON', 0);
addOptional(p, 'maxIter', 500);
addOptional(p, 'rel_err_tol', 1e-6);
addOptional(p, 'loss_tol', 0.0014);
addOptional(p, 'X_true', 0);
addOptional(p, 'operator_name', '');
addOptional(p, 'NesterovON', 1);
addOptional(p, 'U0', 0); 
addOptional(p, 'V0', 0); 
parse(p, A, all_b, r, varargin{:});

maxIter     = p.Results.maxIter;
rel_err_tol = p.Results.rel_err_tol;
loss_tol    = p.Results.loss_tol;
plotON      = p.Results.plotON;
debugON     = p.Results.debugON;
X_true      = p.Results.X_true;
NesterovON  = p.Results.NesterovON;
operator_name         = p.Results.operator_name;
U0 = p.Results.U0;
V0 = p.Results.V0;

[N, ~, M] = size(A);

% lambda0 = 1;
%% some test

% if debugON
%     [U, S, V] = svd(X_true);
%     U_true = U(:, 1:r) * S(1:r, 1:r);
%     V_true = V(:, 1:r);
%     norm(U_true * V_true' - X_true);
% end

if isscalar(U0)
    U0 = randn(N, r) + randn(N, r)*1i;
    V0 = randn(N^2, r) + randn(N^2, r)*1i;
end



%% ALS

ALS_one_loop_counter = tic;
% Initialize U and V randomly

% V0 = V_true;
X0 = U0*V0';

all_U = zeros(N, r, 2);
all_V = zeros(N^2, r, 2);
all_X = zeros(N, N^2, 2);
loss  = zeros(2, 1);

flag = 0;
% ALS Optimization Loop
for iter = 1:maxIter
    % iter
    % iter
    % Step 1: Update U with fixed V
    vec_U1_all = zeros(N*r, N);
    for k = 1:N
        V0_k = V0((k-1)*N+1:k*N, :);
        b_k = all_b(k, :).';

        regmat_U = zeros(M, N*r);
        for m = 1:M
            regmat_U(m, :) = vec(A(:, :, m) * V0_k)'; % Matrix product A_i * V
        end
        vec_U1_all(:, k) = lsqminnorm(regmat_U, b_k);
    end
    
    vec_U1 = mean(vec_U1_all, 2);
    U1_temp = reshape(vec_U1, [N, r]);

    % Step 2: Update V with fixed U
    regmat_V = zeros(M, N*r);
    for m = 1:M
        regmat_V(m, :) = vec(A(:, :, m)' * U1_temp)'; % Matrix product A_i * V
    end

    for k = 1:N
        b_k = all_b(k, :).';
        vec_V1 = lsqminnorm(regmat_V, conj(b_k));
        % vec_V1 = regmat_V\conj(b);
        V1_temp((k-1)*N+1:k*N, :) = reshape(vec_V1, [N, r]);
    end


    % Nesterov-ALS with restarting
    if NesterovON && iter > 3
        beta = 0.4;
        eta = 1.3;

        % Restarting
        if loss(iter-1) > eta * loss(iter-2)
            if debugON
                disp('restart')
            end
            % lambda0 = 1;
            beta = 0;
        end
        % Nesterov momentum updating
        U1 = U1_temp + beta * (all_U(:, :, iter-1) - all_U(:, :, iter-2));
        V1 = V1_temp + beta * (all_V(:, :, iter-1) - all_V(:, :, iter-2));

    else
        U1 = U1_temp;
        V1 = V1_temp;
    end


    % Convergence Check (Frobenius norm)
    X1 = U1_temp * V1_temp'; % Current estimate of M

    all_U(:, :, iter) = U1;
    all_V(:, :, iter) = V1;
    all_X(:, :, iter) = X1;

    b_est = zeros(N, M);
    for k = 1:N
        b_est(k, :) = squeeze(sum(X1(:, (k-1)*N+1:k*N) .* conj(A), [1, 2]));
    end
    loss(iter) = norm(b_est - all_b, 'fro');
    % loss(iter) = norm(squeeze(sum(X1 .*conj(A), [1, 2])) - b);

    rel_error = norm(X1 - X0, 'fro');
    % rel_error = norm(loss(iter) - loss(iter), 'fro');
    % rel_error = loss(iter);
    % if iter > 5 && (rel_error < rel_err_tol || loss(iter) < loss_tol)
    if iter > 5 && (rel_error < rel_err_tol && loss(iter) < loss_tol)
        fprintf('Converged in   %d \titerations with relative error %.12f, loss %.12f\n', iter, rel_error, loss(iter));
        flag = 1;
        break;
    end


    X0 = X1;
    V0 = V1;
    U0 = U1;

    fprintf('iter = %d, loss=%.8f \n', iter, loss(iter))
    
end

M_est = X1;

if flag == 0
    fprintf('Not Converged in   %d \t iterations with relative error %.12f, loss %.12f\n', iter, rel_error, loss(iter));
end
%%
outputInfo.all_U = all_U;
outputInfo.all_V = all_V;
outputInfo.all_X = all_X;
outputInfo.loss = loss;
outputInfo.time = toc(ALS_one_loop_counter);
outputInfo.convergence_flag = flag;
if debugON
    all_error = sqrt(squeeze(sum((abs(all_X - X_true)).^2, [1, 2])));
    figure;hold on;
    plot(log10(loss));
    plot(log10(all_error));
    ttl = [operator_name, ' Final error = ', num2str(all_error(end))];
    title(ttl);
    legend('loss', 'error')

    outputInfo.all_error = all_error;
end




