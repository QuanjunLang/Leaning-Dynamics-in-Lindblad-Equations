function H_part = get_H_part(H_true)

n = length(H_true);

H_part = -1i*(kron(eye(n), H_true) - kron(H_true.', eye(n)));


end