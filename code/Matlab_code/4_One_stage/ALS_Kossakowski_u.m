function u_est = ALS_Kossakowski_u(all_rho_pair_data, v, H, G, p, lambda)



all_rho = squeeze(all_rho_pair_data(:, :, 1, :));
all_drho = squeeze(all_rho_pair_data(:, :, 2, :));
[~, n, M] = size(all_rho);
vt = v';



H_part = get_H_part(H);
A = zeros(M*n^2, (n^2-1)*p);
b = zeros(M*n^2, 1);

for m = 1:M
    A_u = zeros(n^2, (n^2-1)*p);
    rho = all_rho(:, :, m);
    s = 1;
    for q = 1:p
        for k = 1:n^2-1
            temp = sparse(zeros(n^2, 1));
            for l = 1:n^2 - 1
                temp = temp + vt(q, l)*G{k, l}*vec(rho);
            end
            A_u(:, s) = temp;
            s = s+1;
        end
    end

    drho = all_drho(:, :, m);
    b_u = vec(drho) - H_part*vec(rho);
    

    A((m-1)*n^2+1:m*n^2, :) = A_u;
    b((m-1)*n^2+1:m*n^2) = b_u;
end




% constraint = 1;
% lambda = 1;
if ~lambda
    % without positive definite constraint
    vec_u_est = A\b;
    u_est = reshape(vec_u_est, size(v));


else
    % with constraint: (u and v have to be the same)
    bb = A*vec(v) - b;
    x = ridge(bb, A, lambda, 0);
    u_est = reshape(x(2:end) + vec(v), size(v));
end


end