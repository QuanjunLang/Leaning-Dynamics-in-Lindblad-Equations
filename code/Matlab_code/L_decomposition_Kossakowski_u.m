function u_est = L_decomposition_Kossakowski_u(L_vec, v, H, G, n, p, lambda)

A_u = zeros(n^4, (n^2-1)*p);
vt = v';

s = 1;
for q = 1:p
    for k = 1:n^2-1

        temp = sparse(zeros(n^4, 1));
        for l = 1:n^2 - 1
            temp = temp + vt(q, l)*vec(G{k, l});
        end
        A_u(:, s) = temp;s = s+1;
    end
end


H_part_vec = vec(get_H_part(H));
b_u = (L_vec - H_part_vec);

if ~exist('lambda', 'var')
    lambda = 0;
end


% constraint = 1;
% lambda = 1;


if ~lambda
    % without positive definite constraint
    vec_u_est = A_u\b_u;
    u_est = reshape(vec_u_est, size(v));


else
    % with constraint: (u and v have to be the same)
    b = A_u*vec(v) - b_u;
    x = ridge(b, A_u, lambda, 0);
    u_est = reshape(x(2:end) + vec(v), size(v));
end


end