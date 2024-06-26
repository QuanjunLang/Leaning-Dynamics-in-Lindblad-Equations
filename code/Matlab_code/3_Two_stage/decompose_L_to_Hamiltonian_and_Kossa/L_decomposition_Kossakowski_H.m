function [H_est, c_est] = L_decomposition_Kossakowski_H(L, u, v, G, F, n)




K_est = u*v';

% % sanity check
% norm(K - K_est)
% norm(u - u_est)
% norm(v - v_est)


K_part = zeros(size(L));
for k = 1:n^2 - 1
    for l = 1:n^2 - 1
        K_part = K_part + K_est(k, l)*G{k, l};
    end
end






rho = randn(n, n);
L_rho = L*vec(rho);

A_temp = kron(rho.', eye(n, n)) - kron(eye(n, n), rho);

A_temp_proj = zeros(n^2, n^2-1);
for j = 1:n^2 - 1
    A_temp_proj(:, j) = A_temp*vec(F{j});
end

A_R = real(A_temp_proj);
A_I = imag(A_temp_proj);
A = [A_R; A_I];


Lindbladian = K_part*vec(rho);
b_temp = 1i*(L_rho - Lindbladian);

b_temp_reshaped = reshape(b_temp, [], 1);

b_R = real(b_temp_reshaped);
b_I = imag(b_temp_reshaped);
b = [b_R; b_I];

c_est = A\b;
H_est = zeros(n, n);
for i = 1:n^2-1
    H_est = H_est + c_est(i)*F{i};
end



end