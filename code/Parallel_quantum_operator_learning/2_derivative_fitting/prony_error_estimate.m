clc
close all
clear all
addPaths
%%
N = 20;

dx_obs = 1;
dx = 0.01; 
prony_N = 20;
prony_p = 10;
obs_std = 0;



T_obs = prony_N*dx_obs;
T = max(T_obs, T_obs);

lambda_real = -0.1-rand(N, 1);
lambda_imag = -0.1-rand(N, 1)+1i*randn(N, 1);

% lambda = [lambda_real;lambda_imag;conj(lambda_imag)];
lambda = [lambda_real;lambda_imag];
n = length(lambda);
c = randn(n, 1);

xgrid_obs = (0:dx_obs:prony_N*dx_obs)';
xgrid_all = (0:dx:T)';
xgrid = (0:dx:T_obs)';



f = @(x) 0*x;
for i = 1:n
    f = @(x) f(x) + c(i)*exp(lambda(i)*x);
end

% f = @(x) sum(c.*exp(lambda).^(x'))';

df = @(x) 0*x;
for i = 1:n
    df = @(x) df(x) + lambda(i)*c(i)*exp(lambda(i)*x);
end



f_true = squeeze(f(xgrid_obs));
f_obs = f_true+randn(size(xgrid_obs))*obs_std;

H = hankel(f_obs);
figure;plot(log10(svd(H)),'.')






Prony_I.prony_N = prony_N;
Prony_I.prony_p = prony_p;
Prony_I.obs_dx  = dx_obs;
Prony_I.obs_h_grid = f_obs;
Prony_I.polycoef_method = 'MP';
Prony_I.weight_method = 'LS';
Prony_I.root_normalization = 0;
Prony_I.lambda_augmentation = 0;
Prony_I.drop_0 = 0;

prony_result = prony_method(Prony_I);


%%
figure;
subplot(221);hold on;

plot(xgrid_all, real(f(xgrid_all)), 'LineWidth',5)
plot(xgrid_all, real(prony_result.h(xgrid_all)), 'LineWidth',2)
plot(xgrid_obs, real(f(xgrid_obs)), '.', 'MarkerSize',20)

subplot(222);hold on;
plot(xgrid_all, imag(f(xgrid_all)), 'LineWidth',5)
plot(xgrid_all, imag(prony_result.h(xgrid_all)), 'LineWidth',2)
plot(xgrid_obs, imag(f(xgrid_obs)), '.', 'MarkerSize',20)

subplot(223);hold on;
plot(xgrid_all, real(df(xgrid_all)), 'LineWidth',5)
plot(xgrid_all, real(prony_result.dh(xgrid_all)), 'LineWidth',2)
plot(xgrid_obs, real(df(xgrid_obs)), '.', 'MarkerSize',20)

subplot(224);hold on;
plot(xgrid_all, imag(df(xgrid_all)), 'LineWidth',5)
plot(xgrid_all, imag(prony_result.dh(xgrid_all)), 'LineWidth',2)
plot(xgrid_obs, imag(df(xgrid_obs)), '.', 'MarkerSize',20)




L2_err = sum(abs(f(xgrid) - prony_result.h(xgrid)).^2)*dx;
H0_err = sum(abs(df(xgrid) - prony_result.dh(xgrid)).^2)*dx;
H1_err = L2_err + H0_err;

L2_err_obs = sum(abs(f(xgrid_obs) - prony_result.h(xgrid_obs)).^2)*dx_obs;
H0_err_obs = sum(abs(df(xgrid_obs) - prony_result.dh(xgrid_obs)).^2)*dx_obs;
H1_err_obs = L2_err_obs + H0_err_obs;

%%
figure;hold on;
cx = linspace(-1, 1, 100);
plot(cx, sqrt(1-cx.^2),  'b', 'LineWidth',3)
plot(cx, -sqrt(1-cx.^2), 'b', 'LineWidth',3)
axis equal
plot(exp(lambda), 'r.', 'MarkerSize', 20);
plot(exp(prony_result.lam), '+', 'MarkerSize', 10, 'MarkerFaceColor','black');


fprintf('n = %d, Prony_N = %d\n', n, prony_N);
fprintf('Continuous L2 error = %f, H0 error = %f, H1 error = %f\n', L2_err, H0_err, H1_err)
fprintf('Discrete   L2 error = %f, H0 error = %f, H1 error = %f', L2_err_obs, H0_err_obs, H1_err_obs)