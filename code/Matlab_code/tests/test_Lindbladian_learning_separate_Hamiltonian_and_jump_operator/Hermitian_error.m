function H_err = Hermitian_error(H_est, H_true)


n = length(H_est);
H_est = H_est - H_est(1, 1)*eye(n);
H_true = H_true - H_true(1, 1)*eye(n);
H_err = norm(H_est - H_true, 'fro');

% H_est - H_true
% H_err

end