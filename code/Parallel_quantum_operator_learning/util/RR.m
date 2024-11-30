function E = RR(L)

[N, ~] = size(L);
n = sqrt(N);
assert(floor(n) == n)

E = zeros(N, N);
k = 1;

for j = 1:n
    for i = 1:n
        L_ij = L((i-1)*n+1:i*n, (j-1)*n+1:j*n);
        E(:, k) = vec(L_ij);
        k = k + 1;
        % a = 1;
    end
end









end


