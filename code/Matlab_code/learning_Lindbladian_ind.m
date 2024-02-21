function [C_est, C_est_all] = learning_Lindbladian_ind(sysInfo, all_rho, H, C_other, C_true, ind)

% parser = inputParser;
% addRequired(parser, 'sysInfo');
% addRequired(parser, 'all_rho');
% addRequired(parser, 'H');
% addOptional(parser, 'C_other', {});
% % addOptional(parser, 'use_Lindbladian', 0);
% parse(parser,sysInfo, all_rho, varargin{:});
% 
% C_other = parser.Results.C_other;

n = sysInfo.n;
p = sysInfo.p;
dt = sysInfo.dt;
tgrid = sysInfo.tgrid;
M = sysInfo.M;
L = sysInfo.L;

C_ind_current = rand(sysInfo.n, sysInfo.n);

%% Learning Lindbladian

niter = 50;
C_ind_1_est_all = cell(niter, 1);
C_ind_2_est_all = cell(niter, 1);

C_ind_1_est_all{1}  = C_ind_current;
C_ind_2_est_all{1}  = C_ind_current;
C_ind_2_current     = C_ind_2_est_all{1};

all_b = cell(L-1, M);
for m = 1:M
    for l = 1:L-1
        curr_rho = all_rho(:, :, l, m);
        next_rho = all_rho(:, :, l+1, m);
        Lindbladian = combine_Lindbladian(C_other, curr_rho);
        b_temp = (next_rho - curr_rho)/dt  + 1i*(H*curr_rho - curr_rho*H) - Lindbladian;
        all_b{l, m} = reshape(b_temp, [], 1);
    end
end
bb = cat(1, all_b{:});


for j = 1:niter
    % Firstly for one side
    all_A = cell(L-1, M);
    for m = 1:M
        for l = 1:L-1
            curr_rho = all_rho(:, :, l, m);
            A1 = kron((curr_rho*C_ind_2_current).', eye(n, n));
            A2 = kron(eye(n, n), curr_rho*C_ind_2_current);
            A3 = kron(curr_rho.', C_ind_2_current);
            all_A{l, m} = A1 - 0.5*A2 - 0.5*A3;
        end
    end
    AA_1 = cat(1, all_A{:});
    C_ind_1_est = reshape(AA_1\bb, [n, n]);
    C_ind_1_est_all{j} = C_ind_1_est;
    C_ind_1_current = C_ind_1_est;
    
    % For the other side
    all_A = cell(L-1, M);
    for m = 1:M
        for l = 1:L-1
            curr_rho = all_rho(:, :, l, m);
            A1 = kron(eye(n, n), (C_ind_1_current*curr_rho));
            A2 = kron(C_ind_1_current.', curr_rho);
            A3 = kron((C_ind_1_current*curr_rho).', eye(n, n));
            all_A{l, m} = A1 - 0.5*A2 - 0.5*A3;
        end
    end
    AA_2 = cat(1, all_A{:});
    C_ind_2_est = reshape(AA_2\bb, [n, n]);
    C_ind_2_est_all{j} = C_ind_2_est;
    C_ind_2_current = C_ind_2_est;
    
    % if j == 30
    const = mean(C_ind_1_est_all{j}./C_ind_2_est_all{j}', 'all');
    C_est_all{j} = C_ind_1_est_all{j}/sqrt(const);

    % end

    % C_ind_1_est = C_true{1};
    % C_ind_2_est = C_true{1}';

    % C_ind_1_est = C_est_all{j};
    % C_ind_2_est = C_est_all{j}';


    for m = 1:M
        for l = 1:L-1
            curr_rho = all_rho(:, :, l, m);
            next_rho = all_rho(:, :, l+1, m);
            % LL = C_ind_1_est*curr_rho*C_ind_2_est - 0.5*curr_rho*C_ind_2_est*C_ind_1_est - 0.5*C_ind_2_est*C_ind_1_est*curr_rho;
            C_other_and_current = C_other;
            C_other_and_current{end+1} = C_est_all{j};
            LL = combine_Lindbladian(C_other_and_current, curr_rho);
            % LL = 0;
            loss(m, l) = norm((next_rho - curr_rho)/dt  + 1i*(H*curr_rho - curr_rho*H) - LL, 'fro');
        end
    end
    Loss(j) = mean(loss, 'all');
end

%%
% err_1 = zeros(niter, 1);
% err_2 = zeros(niter, 1);
err_3 = zeros(niter, 1);

for j = 1:niter
    % err_1(j) = norm(C_true{1} - C_ind_1_est_all{j}, 'fro');
    % err_2(j) = norm(C_true{1} - C_ind_2_est_all{j}', 'fro');
    err_3(j) = norm(svd(C_true{ind}) - svd(C_est_all{j}));
end

figure;hold on;
% plot(log10(err_1))
% plot(log10(err_2))
plot(log10(err_3))
plot(log10(Loss))
legend('Lindbladian Unitary error', 'Loss')
%%
C_est = C_est_all{end};
end