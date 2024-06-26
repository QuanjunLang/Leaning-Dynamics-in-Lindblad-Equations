function plot_all_prony_modes(all_rho_prony, sysInfo, trueInfo)


figure;hold on;
for ind_1 = 1:sysInfo.n
    for ind_2 = 1:sysInfo.n
        plot(all_rho_prony{ind_1, ind_2}.lam, '.', 'markersize', 5);
    end
end


eigL = eig(trueInfo.L_true);
xL = max(max(abs(real(eigL))), 1);
yL = max(max(abs(imag(eigL))), 1);
xL = 5;
yL = 5;
plot(eigL, '.', 'markersize', 20);
xlim([-xL, xL]); ylim([-yL, yL])
xline(0)
grid on;
title('Eigenvalues of L')

end