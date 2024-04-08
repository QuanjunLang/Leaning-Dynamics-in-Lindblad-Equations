clc
%% partial fractional decomposition

I = I_poly;


Ts = obs_dx;
r = I.r;
w = I.h;
h0 = h_grid_obs(1);

%% Find the inverse of L[h]
syms x
expr = 0;
lam = log(r)/Ts;

for i = 1:p
    expr = expr + w(i)/(x - lam(i));
end

expr_inv = vpa(1/expr); 
simp = partfrac(expr_inv, 'FactorMode', 'complex');
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
    [N, D] = numden(C(i));
    coef = coeffs(D);
    par = coef(2);
    w_seq(i-1) = N/par;
    lam_seq(i-1) = coef(1)/par;
end

theta_hat = @(x) 0*x;
for i = 1:p-1
    theta_hat = @(x) theta_hat(x) + w_seq(i)./(x + lam_seq(i));
end

figure;
subplot(121);hold on;
theta_hat_true = @(x) x.^(-1/2);
plot(xgrid, real(theta_hat(xgrid)))
plot(xgrid, theta_hat_true(xgrid));

xgrid_log = 10.^(linspace(-5, 5, 500));
subplot(122);hold on;
theta_hat_true = @(x) x.^(-1/2);
plot(log10(xgrid_log), real(theta_hat(xgrid_log)));
plot(log10(xgrid_log), theta_hat_true(xgrid_log));

legend('est Lap[gamma]', 'True Lap[gamma]')

%% Inverse Laplace transform (explicit formula)
f = @(x) 0*x;
for i = 1:p-1
    f = @(x) f(x) + w_seq(i)*exp(-lam_seq(i).*x);
end


gamma = @(x) x.^(-1/2)/(sqrt(pi));    
xgrid = 0:0.01:20;
figure; hold on;
plot(xgrid, real(f(xgrid)), '-.', 'MarkerSize', 10)
plot(xgrid, gamma(xgrid));
legend('est', 'true Lap[gamma]')
ylim([0, 4])
    
    
    
   