L = trueInfo.L_true;
E = trueInfo.E_true;
RL = trueInfo.RL_true;
RE = trueInfo.RE_true;


size(observableInfo.O)
size(observableInfo.rho0)
size(all_rho)




switch sysInfo.observable_option
    case 'Multiple_random_observables'
        ind_m = randi([1, sysInfo.M]);
        ind_o = randi([1, sysInfo.N_o]);
        rho0_O_temp = all_rho(ind_o, 1, ind_m);
        rho1_O_temp = all_rho(ind_o, 2, ind_m);

        O_temp = observableInfo.O(:, :, ind_o, ind_m);
        rho0_temp = observableInfo.rho0(:, :, ind_m);

    case 'Single_random_observable'
        ind_m = randi([1, sysInfo.M]);
        rho0_O_temp = all_rho(1, ind_m);
        rho1_O_temp = all_rho(2, ind_m);

        O_temp = observableInfo.O(:, :, ind_m);
        rho0_temp = observableInfo.rho0(:, :, ind_m);
end


rho0_prime_temp = (rho1_O_temp - rho0_O_temp)/sysInfo.dt;


L_error_vec = norm(vec(rho0_temp)' * L' * vec(O_temp) - rho0_prime_temp);
% L_error_vec_prony = norm(vec(rho0_temp)' * L' * vec(O_temp) - all_rho_prony{ind_o, ind_m}.dh(0));
E_error_vec = norm(vec(rho0_temp)' * E' * vec(O_temp) - rho1_O_temp);

fprintf('Lindbladian error, vectorization: %.18f \n', L_error_vec)
% fprintf('Lindbladian error, vectorization prony: %.18f \n', L_error_vec_prony)
fprintf('Channel error, vectorization: %.18f \n', E_error_vec)


%%
UV = kron(conj(O_temp), rho0_temp);
L_error_R = norm(sum(conj(RL).*(UV), 'all') - rho0_prime_temp);
E_error_R = norm(sum(conj(RE).*(UV), 'all') - rho1_O_temp);

fprintf('Lindbladian error, R: %.18f \n', L_error_R)
fprintf('Channel error, R: %.18f \n', E_error_R)



%%




all_rho0    = observableInfo.rho0;
all_O       = observableInfo.O;
[N, ~, M]   = size(all_O);


A = cell(M, 1);
for m = 1:M
    O = all_O(:, :, m);
    rho0 = all_rho0(:, :, m);
    A{m} = kron(conj(O), rho0);
end

A_mat = zeros(N^2, N^2, M);
for m = 1:M
    A_mat(:, :, m) = A{m};
end

b_E = squeeze(all_rho(sysInfo.channel_dt_rate+1, :))';
r_E = trueInfo.rank_RE_true;
[RE_est, outputInfo_E] = ALS(A_mat, b_E, r_E, 'X_true', RE, 'debugON', 1, 'operator_name', 'E');
% norm(sum(conj(A{M}).*trueInfo.RE_true, 'all') - b_E(M))


b_L = squeeze(all_rho(2, :) - all_rho(1, :))'/sysInfo.dt;
r_L = trueInfo.rank_RL_true;
[RL_est, outputInfo_L] = ALS(A_mat, b_L, r_L, 'X_true', RL, 'debugON', 1, 'operator_name', 'L');
% norm(sum(conj(A{M}).*trueInfo.RL_true, 'all') - b_L(M))

% 
% all_E_temp = zeros(size(outputInfo_L.all_M));
% for i = 1:length(all_E_temp)
%     all_E_temp(:, :, i) = RR(expm(RR(outputInfo_L.all_M(:, :, i))*sysInfo.dt*sysInfo.channel_dt_rate));
% end
% all_error_E_temp = squeeze(sum((abs(all_E_temp - trueInfo.RE_true)).^2, [1, 2]));
% all_error_L = squeeze(sum((abs(outputInfo_L.all_M - trueInfo.RL_true)).^2, [1, 2]));
% all_error_E = squeeze(sum((abs(outputInfo_E.all_M - trueInfo.RE_true)).^2, [1, 2]));
% 
% figure;hold on;
% plot(log10(all_error_E))
% plot(log10(all_error_E_temp))
% plot(log10(all_error_L))
% legend('E error', 'E temp error', 'L error')