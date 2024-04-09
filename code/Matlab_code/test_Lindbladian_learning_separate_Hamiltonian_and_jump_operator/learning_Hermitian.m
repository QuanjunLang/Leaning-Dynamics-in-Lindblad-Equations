function H_est = learning_Hermitian(sysInfo, all_rho, varargin)

parser = inputParser;
addRequired(parser, 'sysInfo');
addRequired(parser, 'all_rho');
addOptional(parser, 'C', {});
% addOptional(parser, 'use_Lindbladian', 0);
parse(parser,sysInfo, all_rho, varargin{:});

C = parser.Results.C;

n = sysInfo.n;
p = sysInfo.p;
dt = sysInfo.dt;
tgrid = sysInfo.tgrid;
M = sysInfo.M;
L = sysInfo.L;

if length(C) > 0
    use_Lindbladian = 1;
else
    use_Lindbladian = 0;
end


%% Learning Hermitian
all_A = cell(L-1, M);
all_b = cell(L-1, M);
for m = 1:M
    for l = 1:L-1
        curr_rho = all_rho(:, :, l, m);
        next_rho = all_rho(:, :, l+1, m);
        all_A{l, m} = kron(curr_rho.', eye(n, n)) - kron(eye(n, n), curr_rho);

        if use_Lindbladian
            Lindbladian = combine_Lindbladian(C, curr_rho);
            b_temp = 1i*((next_rho - curr_rho)/dt - Lindbladian);
        else
            b_temp = (next_rho - curr_rho)/dt * 1i;
        end
        all_b{l, m} = reshape(b_temp, [], 1);
    end
end

AA = cat(1, all_A{:});
bb = cat(1, all_b{:});

%% Trivial identifiability issue
% H + CI is still a valid result
% make the leading coef to be 0

H_est = reshape(AA\bb, [n, n]);






end