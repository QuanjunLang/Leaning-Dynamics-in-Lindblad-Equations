function E = rearrangement_R(L)

[N, ~] = size(L);
n = sqrt(N);
assert(floor(n) == n)

E = zeros(N, N);
k = 1;

for i = 1:n
    for j = 1:n
        L_ij = L((i-1)*n+1:i*n, (j-1)*n+1:j*n);
        E(k, :) = reshape(L_ij.', [], 1);
        k = k + 1;
        % a = 1;
    end
end









end


