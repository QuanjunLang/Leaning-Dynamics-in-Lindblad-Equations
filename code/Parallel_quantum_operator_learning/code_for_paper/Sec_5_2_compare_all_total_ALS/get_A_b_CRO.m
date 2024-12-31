function [A, b] = get_A_b_CRO(sysInfo, observableInfo, RK)


% Fixed first row Observables
n = sysInfo.n;
M = sysInfo.M;
N_o = sysInfo.N_o;
all_rho0 = observableInfo.rho0;


random_seed = randi(10000);
a = pyrunfile("Generate_observables.py", 'a', n = n, N_o = N_o*M, random_seed=random_seed);
many_random_observables = double(a);

A = zeros(n^2, n^2, M*N_o);
b = zeros(M*N_o, 1);

l = 1;
for m = 1:M
    rho0 = all_rho0(:, :, m);
    for j = 1:N_o
        O = squeeze(many_random_observables(l, :, :));
        
        A(:, :, l) = kron(conj(O), rho0);
        b(l) = sum(conj(A(:, :, l)).*RK, 'all');
        l = l+1;
    end
end


end