
L = trueInfo.L_true;
E = trueInfo.E_true;
RL = trueInfo.RL_true;
RE = trueInfo.RE_true;

if ~ sysInfo.UseDiffObs
    ind_m = sysInfo.M;
    ind_o = sysInfo.N_o;

    O_temp = observableInfo.O(:, :, ind_o, ind_m);
    rho0_temp = observableInfo.rho0(:, :, ind_m);

    rho0_O_temp = all_rho(ind_o, 1, ind_m);
    rho1_O_temp = all_rho(ind_o, 2, ind_m);

else
    ind_m = sysInfo.M;

    O_temp = observableInfo.O(:, :, ind_m);
    rho0_temp = observableInfo.rho0(:, :, ind_m);

    rho0_O_temp = all_rho(1, ind_m);
    rho1_O_temp = all_rho(2, ind_m);
end

rho0_prime_temp = (rho1_O_temp - rho0_O_temp)/sysInfo.dt;


L_error_vec = norm(vec(rho0_temp)' * L' * vec(O_temp) - rho0_prime_temp);
% L_error_vec_prony = norm(vec(rho0_temp)' * L' * vec(O_temp) - all_rho_prony{ind_o, ind_m}.dh(0));
E_error_vec = norm(vec(rho0_temp)' * E' * vec(O_temp) - rho1_O_temp);

fprintf('Lindbladian error, vectorization: %.18f \n', L_error_vec)
% fprintf('Lindbladian error, vectorization prony: %.18f \n', L_error_vec_prony)
fprintf('Channel error, vectorization: %.18f \n', E_error_vec)

UV = kron(conj(rho0_temp), O_temp);
L_error_R = norm(sum(conj(RL).*(UV), 'all') - rho0_prime_temp);
E_error_R = norm(sum(conj(RE).*(UV), 'all') - rho1_O_temp);

fprintf('Lindbladian error, R: %.18f \n', L_error_R)
fprintf('Channel error, R: %.18f \n', E_error_R)


L_innerp_sym_error = sum(conj(RL).*(UV), 'all') - sum(conj(UV) .* RL, 'all');
E_innerp_sym_error = sum(conj(RE).*(UV), 'all') - sum(conj(UV) .* RE, 'all');

L_sym = norm(RL - RL');
E_sym = norm(RE - RE');

fprintf('L inner product sym error %.18f \n', L_innerp_sym_error)
fprintf('E inner product sym error %.18f \n', E_innerp_sym_error)

fprintf('L inner product sym error %.18f \n', L_sym)
fprintf('E inner product sym error %.18f \n', E_sym)