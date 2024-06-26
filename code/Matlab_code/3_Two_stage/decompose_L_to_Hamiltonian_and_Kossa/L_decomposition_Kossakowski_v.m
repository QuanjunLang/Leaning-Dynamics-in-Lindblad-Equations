function v_est = L_decomposition_Kossakowski_v(L_vec, u, H, G, n, p, lambda)


A_v = zeros(n^4, (n^2-1)*p);
% u = v_true;
% v = v_true;

s = 1;
for q = 1:p
    for l = 1:n^2-1
        temp = sparse(zeros(n^4, 1));
        for k = 1:n^2 - 1
            temp = temp + u(k, q)*vec(G{k, l});
        end
        A_v(:, s) = temp;s = s+1;
    end
end


H_part_vec = vec(get_H_part(H));
b_v = (L_vec - H_part_vec);



A_v = conj(A_v);
b_v = conj(b_v);

if ~exist('lambda', 'var')
    lambda = 0;
end



% constraint = 0;
% lambda = 1;
if ~lambda
    % without positive definite constraint
    vec_v_est = A_v\b_v;
    v_est = reshape(vec_v_est, size(u));

else
    % with constraint: (u and v have to be the same)
    b = A_v*vec(u) - b_v;
    x = ridge(b, A_v, lambda, 0);
    v_est = reshape(x(2:end) + vec(u), size(u));
end

% % v_est = conj(v_est);





end