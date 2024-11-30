function result = multi_mat_als_observable(expect_derivative_0, observableInfo, r, varargin)
% Use ALS to learn the generator of Lindbladian Master Equation from observables

%% load parameters
p = inputParser;

addRequired(p, 'expect_derivative_0');
addRequired(p, 'r');
addRequired(p, 'observableInfo');

addOptional(p, 'trueInfo', 0);
addOptional(p, 'plotON', 0);
addOptional(p, 'niter', 2000);
addOptional(p, 'loss_threshold', 1e-10);

parse(p, expect_derivative_0, observableInfo, r, varargin{:});

trueInfo           = p.Results.trueInfo;
plotON             = p.Results.plotON;
niter              = p.Results.niter;
loss_threshold     = p.Results.loss_threshold;


all_O       = observableInfo.O;
all_rho_0   = observableInfo.rho0;
[N_o, M]    = size(expect_derivative_0);
[n, ~, ~]      = size(observableInfo.rho0);

%% Initialize
niter_pre = 2;
all_A = zeros(n, n, r, niter_pre);      % All steps for A
all_B = zeros(n, n, r, niter_pre);      % All steps for B
all_L = zeros(n^2, n^2, niter_pre);     % All steps for L
all_E = zeros(n^2, n^2, niter_pre);

loss  = zeros(niter_pre, 1);
all_rk = zeros(niter_pre, 1);

A_0 = randn(n, n, r) + 1i*randn(n, n, r);
B_0 = randn(n, n, r) + 1i*randn(n, n, r);

if isa(trueInfo, 'struct')             % If true para information is provided
    L_true = trueInfo.L_true;
    E_true = trueInfo.E_true;

    err_A = zeros(niter_pre, 1);        % store the errors
    err_B = zeros(niter_pre, 1);
    err_V = zeros(niter_pre, 1);
    err_L = zeros(niter_pre, 1);
    err_E = zeros(niter_pre, 1);
end

% %% Sanity check: Using true parameters as initial values
% [U, S, V] = svd(trueInfo.E_true);
% 
% UU = U(:, 1:r)*S(1:r, 1:r);
% VV = conj(V(:, 1:r));
% 
% BB = (permute(reshape(U(:, 1:r)*S(1:r, 1:r), [n, n, r]), [1, 2, 3]));
% AA = (permute(reshape(conj(V(:, 1:r)), [n, n, r]), [2, 1, 3]));
% 
% EE = zeros(n^2, n^2);
% for k = 1:r
%     EE = EE + vec(BB(:, :, k))*vec(permute(AA(:, :, k), [2, 1])).';
% end
% 
% 
% fprintf('difference btw E true and svd constructed: %d\n', sum((EE - trueInfo.E_true).^2, 'all'))
% 
% A_0 = AA;
% B_0 = BB;
% 
% 
% A = A_0(:, :, 1);
% B = B_0(:, :, 1);
% %% L
% L_temp = zeros(n^2, n^2);
% for i = 1:r
%     L_temp = L_temp + kron(B_0(:, :, i).', A_0(:, :, i));
% end
% 
% fprintf('difference btw L true and svd constructed: %d\n', sum((L_temp - L_true).^2, 'all'))
% 
% 
% %% vec
% rho0  = all_rho_0(:, :, 1);
% O_cur = all_O(:, :, 1, 1);
% rhod0 = expect_derivative_0(1, 1);
% L_rho_0 = zeros(n, n);
% for k = 1:r
%     L_rho_0 = L_rho_0 + A_0(:, :, k) * rho0 * B_0(:, :, k);
% end
% 
% fprintf('difference btw vec first and vec last: %d\n', sum((vec(L_rho_0) - L_true*vec(rho0)).^2, 'all'))
% 
% 
% %% Observable
% fprintf('difference btw observable (mat:trace): %d\n', abs(trace(O_cur' * L_rho_0) - rhod0))
% fprintf('difference btw observable (mat:Frobn): %d\n', abs(sum(O_cur.' .* L_rho_0, 'all') - rhod0))
% 
% 
% rhod0_vec = vec(rhod0);
% O_cur_vec = vec(O_cur.');% There has to be a .'
% rho0_vec = vec(rho0);
% 
% 
% fprintf('difference btw observable (vec): %d\n', abs((L_true*rho0_vec).' * O_cur_vec - rhod0_vec))
% 
% 
% 
% %% estimate B one line
% % (L_true*vec(rho0)).'* O_cur_vec
% 
% B_temp = 0;
% for k = 1:r
%     B_temp = B_temp + O_cur_vec.' *kron(eye(n, n), A_0(:, :, k)*rho0) * vec(BB(:, :, k));
% end
% 
% % B_temp
% 
% 
% 
% 





%% Iterations

for i = 1:niter
    % Estimate B first
    % tic
    % M = 1;
    regmat_B_given_A = zeros(M*N_o, r*n^2);
    % bB = zeros(M*N_o, 1);
    for m = 1:M
        rho0 = all_rho_0(:, :, m);
        PP = superkron(eye(n, n), pagemtimes(A_0, rho0));
        PPA = reshape(PP, [n^2, n^2*r]);

        O_mat = reshape(all_O(:, :, :, m), [n^2, N_o])';


        regmat_B_given_A((m-1)*N_o+1:m*N_o, :) = O_mat*PPA;
        % bB((m-1)*N_o+1:m*N_o) = expect_derivative_0(:, m);
    end
    b = vec(expect_derivative_0);

    vec_B_est = regmat_B_given_A \ b;
    B_0 = reshape(vec_B_est, [n, n, r]);

    % B_true = reshape(BB, [], 1);
    % regmat_B_given_A*B_true - bB




    all_B(:, :, :, i) = B_0;

    % Estimate A next
    regmat_A_given_B = zeros(M*N_o, r*n^2);
    for m = 1:M
        rho0 = all_rho_0(:, :, m);
        % PP = superkron(eye(n, n), pagemtimes(permute(B_0, [2, 1, 3]), permute(rho0, [2,1,3,4])));
        PP = superkron(permute(pagemtimes(rho0, B_0), [2, 1, 3]), eye(n, n));
        PPA = reshape(PP, [n^2, n^2*r]);

        O_mat = reshape(all_O(:, :, :, m), [n^2, N_o])';


        regmat_A_given_B((m-1)*N_o+1:m*N_o, :) = O_mat*PPA;
    end

    vec_A_est = regmat_A_given_B \ b;

    A_0 = reshape(vec_A_est, [n, n, r]);
    all_A(:, :, :, i) = A_0;






    % % balance A and B to contruct V
    % Const = A_0 ./ conj(permute(B_0, [2, 1, 3]));
    % V_0 = A_0./sqrt(Const);
    % all_V(:, :, :, i) = V_0;

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
    % all_E_svd(:, i) = svd(E_0);


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
    if isa(trueInfo, 'struct')
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
    expect_derivative_0_est = zeros(N_o, M);
    for m = 1:M
        rho0 = all_rho_0(:, :, m);
        temp = pagemtimes(A_0, rho0);
        drho = sum(pagemtimes(temp, B_0), 3);
        expect_derivative_0_est(:, m) = squeeze(sum(drho.'.*all_O(:, :, :, m), [1, 2]));
    end
    loss(i) = norm(expect_derivative_0_est - expect_derivative_0, 'fro');


    






    if i > 2
        if abs(loss(i)-loss(i-1)) < loss_threshold
            break
        end
    end

    if mod(i, 100)== 1
        fprintf('niter = %d, loss = %f \n', i, loss(i))
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
% result.all_V = all_V;
result.all_E = all_E;
result.all_rk = all_rk;

if isa(trueInfo, 'struct')
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
    if isa(trueInfo, 'struct')
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