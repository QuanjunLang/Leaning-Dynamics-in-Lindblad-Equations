function result = multi_mat_als(all_rho, true_para, r)

[n, ~, T_temp, M] = size(all_rho);
T = T_temp - 1;

% r = 5;
niter = 20;

all_A = zeros(n, n, r, niter);      % All steps for A
all_B = zeros(n, n, r, niter);      % All steps for B
all_V = zeros(n, n, r, niter);      % All steps for V, V is a formal average of A and B', because that is what they are
all_K = zeros(n^2, n^2, niter);     % All steps for K
all_E = zeros(n^2, n^2, niter);

err_A = zeros(niter, 1);        % store the errors
err_B = zeros(niter, 1);
err_V = zeros(niter, 1);
err_K = zeros(niter, 1);
err_E = zeros(niter, 1);
loss  = zeros(niter, 1);

A_0 = randn(n, n, r);
B_0 = randn(n, n, r);

A_true = true_para.A_true;
B_true = true_para.B_true;
r_true = true_para.r_true;
K_true = true_para.K_true;
E_true = true_para.E_true;


% A_0 = A_true; % sanity check of starting from truth
for i = 1:niter
    % Estimate B first
    regmat_B_given_A = zeros(M*T*n^2, r*n^2);
    k = 1;
    for m = 1:M
        for t = 1:T
            rho0 = all_rho(:, :, t, m);
            PP = superkron(eye(n, n), pagemtimes(A_0, rho0));
            PPA = reshape(PP, [n^2. n^2*r]);
            regmat_B_given_A((k-1)*n^2+1:k*n^2, :) = PPA;
            k = k+1;
        end
    end
    bB = vec(all_rho(:, :, 2:end, :));
    vec_B_est = regmat_B_given_A \ bB;
    B_0 = reshape(vec_B_est, [n, n, r]);
    all_B(:, :, :, i) = B_0;

    % Estimate A next
    regmat_A_given_B = zeros(M*T*n^2, r*n^2);
    k = 1;
    for m = 1:M
        for t = 1:T
            rho0 = all_rho(:, :, t, m);
            PP = superkron(eye(n, n), pagemtimes(permute(B_0, [2, 1, 3]), permute(rho0, [2,1,3,4])));
            PPA = reshape(PP, [n^2. n^2*r]);
            regmat_A_given_B((k-1)*n^2+1:k*n^2, :) = PPA;
            k = k+1;
        end
    end

    bA = vec(permute(all_rho(:, :, 2:end, :), [2,1,3,4]));
    vec_A_est = regmat_A_given_B \ bA;
    A_0 = permute(reshape(vec_A_est, [n, n, r]), [2,1,3]);
    all_A(:, :, :, i) = A_0;

    % balance A and B to contruct V
    Const = A_0 ./ conj(permute(B_0, [2, 1, 3]));
    V_0 = A_0./sqrt(Const);
    all_V(:, :, :, i) = V_0;
    
    % contruct K
    K_0 = zeros(n^2, n^2);
    for k = 1:r
        K_0 = K_0 + kron(B_0(:, :, k).', A_0(:, :, k));
    end
    all_K(:, :, i) = K_0;

    % contruct E
    E_0 = zeros(n^2, n^2);
    for k = 1:r
        E_0 = E_0 + vec(B_0(:, :, k))*vec(permute(A_0(:, :, k), [2, 1])).';
    end
    all_E(:, :, i) = E_0;
    
    % Compute error at each step
    err_K(i) = norm(K_0 - K_true, 'fro');
    err_E(i) = norm(E_0 - E_true, 'fro');
    % Note that it is meaningless to compare A_true and A_0 with different r

    if r == r_true
        err_A(i) = norm(A_0 - A_true, 'fro');
        err_B(i) = norm(B_0 - B_true, 'fro');
        % Error of V has to be computed in the sense of SVD because of the
        % unitary invarnce, this works when r = 1
        % When r is large, this is not enough, since V are not unique
        err_V(i) = 0;
        for k = 1:r
            err_V(i) = err_V(i) + norm(svd(V_0(:, :, k)) - svd(A_true(:, :, k)), 'fro');
        end
    end


    % This is the loss, showing that we are really minimizing the loss
    test_rho = zeros(size(all_rho));
    test_rho(:, :, 1, :) = all_rho(:, :, 1, :);
    for t = 1:T
        all_rho0 = all_rho(:, :, t, :);
        all_rho1 = zeros(n, n, 1, M);
        for k = 1:r
            temp = pagemtimes(A_0(:, :, k), all_rho0);
            all_rho1 = all_rho1 + pagemtimes(temp, B_0(:, :, k));
        end
        test_rho(:, :, t+1, :) = all_rho1;
    end
    loss(i) = norm(test_rho - all_rho, 'fro');
end


%% decompose E_est to mat A and B
E_est = all_E(:, :, end);
E_est = (E_est + E_est')/2;
r_est = rank(E_est);
[U, ~, V] = svd(E_est);

vec_V = (U(:, 1:3) + V(:, 1:3))/2;
V_est = zeros(n, n, r_est);
for k = 1:r_est
    V_est(:, :, k) = reshape(vec_V(:, k), [n, n]);
end

result.all_A = all_A;
result.all_B = all_B;
result.all_K = all_K;
result.all_V = all_V;
result.all_E = all_E;

result.err_A = err_A;
result.err_B = err_B;
result.err_V = err_V;
result.err_K = err_K;
result.err_E = err_E;
result.loss  = loss;

result.K_est = all_K(:, :, end);
result.E_est = E_est;
result.V_est = V_est;

%%
figure;hold on;
plot(log10(loss), 'DisplayName','Loss');
plot(log10(err_V), 'DisplayName','V');
plot(log10(err_B), 'DisplayName','B');
plot(log10(err_A), 'DisplayName','A');
plot(log10(err_K), '-', 'linewidth', 3, 'DisplayName','K');
plot(log10(err_E), ':', 'linewidth', 3, 'DisplayName','E');
legend()


end