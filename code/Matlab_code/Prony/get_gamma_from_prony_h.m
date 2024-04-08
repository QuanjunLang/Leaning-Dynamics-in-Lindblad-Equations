function result_gamma = get_gamma_from_prony_h(I, result_prony)
%% partial fractional decomposition

% Ts = I.obs_dx;
% xgrid = I.xgrid;

% r = result_prony.r;
w = result_prony.w;
lam = result_prony.lam;
h0 = result_prony.h(0);


p = I.prony_p;
%% Find the inverse of L[h]
syms x
expr = 0;

for i = 1:p
    expr = expr + w(i)/(x - lam(i));
end

expr_inv = vpa(1/expr);
simp = h0*partfrac(expr_inv, 'FactorMode', 'complex');
% pretty(vpa(simp))

%% Get the laplace transform of memory kernel
C = children(simp);
% for i = 1:length(C)
%     vpa(C(i))
% end

temp = 0;
w_seq = zeros(length(C) - 2, 1);
lam_seq = zeros(length(C) - 2, 1);

for i = 2:length(C) - 1
    [N, D] = numden(C{i});
    coef = coeffs(D);
    par = coef(2);
    w_seq(i-1) = N/par;
    lam_seq(i-1) = coef(1)/par;
end
lam_seq = -lam_seq;
const = C(end);
%% singularity of gamma ------ negative lambda
for i = 1:length(lam_seq)
    if real(lam_seq(i)) > 0
        lam_seq(i) = -1e-5 + imag(lam_seq(i))*1i;
    end
end


%%
theta_hat = @(x) 0*x + const;
for i = 1:length(w_seq)
    theta_hat = @(x) theta_hat(x) + w_seq(i)./(x - lam_seq(i));
end

%% Inverse Laplace transform (explicit formula)
f = @(x) 0*x;
% f = @(x) const*normpdf(x, 0, 1e-3);
for i = 1:length(w_seq)
    f = @(x) f(x) + w_seq(i)*exp(lam_seq(i).*x);
end

result_gamma.lam_seq = lam_seq;
result_gamma.w_seq = w_seq;
result_gamma.gamma = f;
result_gamma.gamma_hat = theta_hat;

S = @(x) 0*x;
for k = 1:length(w_seq)
    S = @(x) S(x) + w_seq(k)*(-lam_seq(k))/pi ./(x.^2 + lam_seq(k)^2);
end

result_gamma.gamma_fourier = S;

plotON = 0;
if plotON
    
    
    
    
    figure;
    subplot(121);hold on;
    theta_hat_true = @(x) x.^(-1/2);
    plot(xgrid, real(theta_hat(xgrid)), 'DisplayName', 'est Lap[gamma]');
    plot(xgrid, theta_hat_true(xgrid), 'DisplayName', 'True Lap[gamma]');
    legend;
    
    subplot(122);hold on;
    theta_hat_true = @(x) x.^(-1/2);
    xgrid_log = 10.^(linspace(-5, 5, 500));
    plot(log10(xgrid_log), real(theta_hat(xgrid_log)), 'DisplayName', 'est Lap[gamma]');
    plot(log10(xgrid_log), theta_hat_true(xgrid_log), 'DisplayName', 'True Lap[gamma]');
    legend;
    
    
    figure;
    subplot(121);hold on;
    
    gamma = @(x) x.^(-1/2)/(sqrt(pi));
    xgrid = 0:0.01:20;
    
    plot(xgrid, real(f(xgrid)), '-.', 'MarkerSize', 10, 'displayname', 'est')
    plot(xgrid, gamma(xgrid), 'displayname', 'true Lap[gamma]')
    legend;
    ylim([0, 4])
    title('Memory kernel')
    
    subplot(122);hold on;grid on;
    plot(log10(xgrid_log), log10(real(f(xgrid_log))), '-.', 'MarkerSize', 10, 'displayname', 'est')
    plot(log10(xgrid_log), log10(gamma(xgrid_log)), 'displayname', 'true Lap[gamma]')
    legend;
    ylim([-2, 4])
    title('Memory kernel, log scale')
    
end
end
