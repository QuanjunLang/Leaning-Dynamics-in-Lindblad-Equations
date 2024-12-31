function [A_FFO, b_FFO] = get_A_b_FFO(sysInfo, observableInfo, all_rho)


% Fixed first row Observables
n = sysInfo.n;
M = sysInfo.M;
N_o = sysInfo.N_o;
all_rho0 = observableInfo.rho0;


A_FFO = zeros(n^2, n^2, M*N_o);
l = 1;
for m = 1:M
    rho0 = all_rho0(:, :, m);
    for j = 1:N_o
        O = observableInfo.O(:, :, j);
        A_FFO(:, :, l) = kron(conj(O), rho0);
        l = l+1;
    end
end

b_FFO = vec(all_rho);


end