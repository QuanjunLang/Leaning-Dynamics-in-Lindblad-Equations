function [H_est, c_est] = ALS_Kossakowski_H(all_rho_pair_data, u, v, G, F)

all_rho = squeeze(all_rho_pair_data(:, :, 1, :));
all_drho = squeeze(all_rho_pair_data(:, :, 2, :));
[~, n, M] = size(all_rho);


K_est = u*v';

% % sanity check
% norm(K - K_est)
% norm(u - u_est)
% norm(v - v_est)


K_part = zeros(n^2, n^2);
for k = 1:n^2 - 1
    for l = 1:n^2 - 1
        K_part = K_part + K_est(k, l)*G{k, l};
    end
end

A = zeros(2*M*n^2, n^2-1);
b = zeros(2*M*n^2, 1);

for m = 1:M
    rho = all_rho(:, :, m);
    drho = all_drho(:, :, m);

    A_temp = kron(rho.', eye(n, n)) - kron(eye(n, n), rho);

    A_temp_proj = zeros(n^2, n^2-1);
    for j = 1:n^2 - 1
        A_temp_proj(:, j) = A_temp*vec(F{j});
    end

    A_R = real(A_temp_proj);
    A_I = imag(A_temp_proj);
    A_m = [A_R; A_I];


    Lindbladian = K_part*vec(rho);
    b_temp = 1i*(vec(drho) - Lindbladian);

    b_temp_reshaped = reshape(b_temp, [], 1);

    b_R = real(b_temp_reshaped);
    b_I = imag(b_temp_reshaped);
    b_m = [b_R; b_I];
    
    A((m-1)*2*n^2+1:m*n^2*2, :) = A_m;
    b((m-1)*2*n^2+1:m*n^2*2) = b_m;
end





c_est = A\b;
H_est = zeros(n, n);
for i = 1:n^2-1
    H_est = H_est + c_est(i)*F{i};
end



end