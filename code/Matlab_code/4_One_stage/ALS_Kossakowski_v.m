function v_est = ALS_Kossakowski_v(all_rho_pair_data, u, H, G, p, lambda)



all_rho = squeeze(all_rho_pair_data(:, :, 1, :));
all_drho = squeeze(all_rho_pair_data(:, :, 2, :));
[~, n, M] = size(all_rho);




H_part = get_H_part(H);
A = zeros(M*n^2, (n^2-1)*p);
b = zeros(M*n^2, 1);
for m = 1:M
    A_v = zeros(n^2, (n^2-1)*p);
    rho = all_rho(:, :, m);
    s = 1;
    for q = 1:p
        for l = 1:n^2-1
            temp = sparse(zeros(n^2, 1));
            for k = 1:n^2 - 1
                temp = temp + u(k, q)*G{k, l}*vec(rho);
            end
            A_v(:, s) = temp;
            s = s+1;
        end
    end


    drho = all_drho(:, :, m);
    b_v = vec(drho) - H_part*vec(rho);
    

    A((m-1)*n^2+1:m*n^2, :) = A_v;
    b((m-1)*n^2+1:m*n^2) = b_v;
end
% 
% H_part_vec = vec(get_H_part(H));
% b_v = (L_vec - H_part_vec);



A = conj(A);
b = conj(b);



% constraint = 0;
% lambda = 1;
if ~lambda
    % without positive definite constraint
    vec_v_est = A\b;
    v_est = reshape(vec_v_est, size(u));

else
    % with constraint: (u and v have to be the same)
    bb = A*vec(u) - b;
    x = ridge(bb, A, lambda, 0);
    v_est = reshape(x(2:end) + vec(u), size(u));
end

% % v_est = conj(v_est);





end