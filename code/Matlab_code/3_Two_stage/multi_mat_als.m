function result = multi_mat_als(all_rho, r, varargin)


%%
p = inputParser;

addRequired(p, 'all_rho');
addRequired(p, 'r');
addOptional(p, 'true_para', 0); 
addOptional(p, 'rank_threshold', 1e-6);
addOptional(p, 'plotON', 0);

parse(p,all_rho, r, varargin{:});

true_para          = p.Results.true_para;
rank_threshold     = p.Results.rank_threshold;
plotON             = p.Results.plotON;



loss_threshold = 1e-10;


[n, ~, T_temp, M] = size(all_rho);
T = T_temp - 1;

% r = 5;
niter_pre = 2;

all_A = zeros(n, n, r, niter_pre);      % All steps for A
all_B = zeros(n, n, r, niter_pre);      % All steps for B
all_V = zeros(n, n, r, niter_pre);      % All steps for V, V is a formal average of A and B', because that is what they are
all_L = zeros(n^2, n^2, niter_pre);     % All steps for L
all_E = zeros(n^2, n^2, niter_pre);

% all_A = cell(niter, 1);
% all_B = cell(niter, 1);
% all_V = cell(niter, 1);



loss  = zeros(niter_pre, 1);
all_rk = zeros(niter_pre, 1);


A_0 = randn(n, n, r) + 1i*randn(n, n, r);
B_0 = randn(n, n, r) + 1i*randn(n, n, r);


if isa(true_para, 'struct')              % If true para information is provided
    % A_true = true_para.A_true;
    % B_true = true_para.B_true;
    r_true = true_para.r_true;
    L_true = true_para.L_true;
    E_true = true_para.E_true;

    err_A = zeros(niter_pre, 1);        % store the errors
    err_B = zeros(niter_pre, 1);
    err_V = zeros(niter_pre, 1);
    err_L = zeros(niter_pre, 1);
    err_E = zeros(niter_pre, 1);

end



niter = 500;
% A_0 = A_true; % sanity check of starting from truth
for i = 1:niter
    i
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

    % contruct L
    L_0 = zeros(n^2, n^2);
    for k = 1:r
        L_0 = L_0 + kron(B_0(:, :, k).', A_0(:, :, k));
    end
    all_L(:, :, i) = L_0;
    
    % contruct E
    E_0 = zeros(n^2, n^2);
    for k = 1:r
        E_0 = E_0 + vec(B_0(:, :, k))*vec(permute(A_0(:, :, k), [2, 1])).';
    end
    all_E(:, :, i)  = E_0;
    all_rk(i)       = rank(E_0);
    all_E_svd(:, i) = svd(E_0);
    
    
    % %% reduce the rank 
    % [UU, S, VV] = svd(E_0);
    % r_updated = sum(diag(S)>rank_threshold);
    % if r_updated < r
    %     BB_0 = UU(:, 1:r_updated)*S(1:r_updated, 1:r_updated);
    %     AA_0 = VV(:, 1:r_updated);
    % 
    %     A_0 = permute(reshape(AA_0, [n, n, r_updated]), [2, 1, 3]);
    %     B_0 = reshape(BB_0, [n, n, r_updated]);
    % end


    %% Compute error at each step, if true parameters are provided
    if isa(true_para, 'struct')
        err_L(i) = norm(L_0 - L_true, 'fro')/n;
        err_E(i) = norm(E_0 - E_true, 'fro')/n;
        % Note that it is meaningless to compare A_true and A_0 with different r

        % if r == r_true
        %     err_A(i) = norm(A_0 - A_true, 'fro');
        %     err_B(i) = norm(B_0 - B_true, 'fro');
        %     % Error of V has to be computed in the sense of SVD because of the
        %     % unitary invarnce, this works when r = 1
        %     % When r is large, this is not enough, since V are not unique
        %     err_V(i) = 0;
        %     for k = 1:r
        %         err_V(i) = err_V(i) + norm(svd(V_0(:, :, k)) - svd(A_true(:, :, k)), 'fro');
        %     end
        % end
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

    if i > 2
        if abs(loss(i) - loss(i-1)) < loss_threshold
            break
        end
    end
            
end


%% decompose E_est to mat A and B
E_est = all_E(:, :, end);
% r_est = all_rk(end);

% If E is symmetric

% E_est = (E_est + E_est')/2;
r_est = rank(E_est);


[U, S, V] = svd(E_est);

vec_V = (U(:, 1:r_est) + V(:, 1:r_est))/2;
V_est = zeros(n, n, r_est);
for k = 1:r_est
    V_est(:, :, k) = reshape(vec_V(:, k), [n, n]);
end

result.all_A = all_A;
result.all_B = all_B;
result.all_L = all_L;
result.all_V = all_V;
result.all_E = all_E;
result.all_rk = all_rk;

if isa(true_para, 'struct') 
    result.err_A = err_A;
    result.err_B = err_B;
    result.err_V = err_V;
    result.err_L = err_L;
    result.err_E = err_E;
end

result.loss  = loss;

result.L_est = all_L(:, :, end);
result.E_est = E_est;
result.V_est = V_est;
result.S_est = diag(S);
%%
if plotON
    figure;
    
    subplot(221)
    hold on;
    plot(log10(loss), 'DisplayName','Loss');
    if isa(true_para, 'struct')
        plot(log10(err_V), 'DisplayName','V');
        plot(log10(err_B), 'DisplayName','B');
        plot(log10(err_A), 'DisplayName','A');
        plot(log10(err_L), '-', 'linewidth', 3, 'DisplayName','L');
        plot(log10(err_E), ':', 'linewidth', 3, 'DisplayName','E');
    end
    legend()
    
    subplot(222)
    plot(log10(all_E_svd)')

end

end