
% close all



w = result_prony.w;
r = result_prony.r;
lam = result_prony.lam;
p = length(w);


dx = 0.01;
L = 10;
xgrid = 0:dx:L;
N = length(xgrid);

%%
h = @(x) real(result_prony.h(x));

g = @(x) 0*x;
for t = 1:p
    g = @(x) g(x) + lam(t).*w(t).*exp(lam(t)*x);
end
g = @(x) real(g(x));

dg = @(x) 0*x;
for t = 1:p
    dg = @(x) dg(x) + lam(t).*lam(t).*w(t).*exp(lam(t)*x);
end

% %%
% I.data_type = 'grid';
% I = update_I(I);
% 
% %%
% figure;
% subplot(131);hold on;
% plot(I.xgrid, prony_h_grid);
% plot(I.xgrid, I.h(I.xgrid));
% plot(I.xgrid, I.h_grid);
% title("h")
% legend('True', 'Prony', 'noisy')
% 
% %%
% subplot(132);hold on;
% plot(I.xgrid, gradient(I.h(I.xgrid), dx));
% plot(I.xgrid, g(I.xgrid));
% plot(I.xgrid, gradient(I.h_grid, dx));
% title("h'")
% legend('True', 'Prony', 'noisy')
% 
% %%
% subplot(133);hold on;
% plot(I.xgrid, gradient(gradient(I.h(I.xgrid), dx), dx));
% plot(xgrid, dg(xgrid));
% plot(I.xgrid, gradient(gradient(I.h_grid, dx), dx));
% title("h''")
% legend('True', 'Prony', 'noisy')
%%
K = zeros(N, N);
for i = 1:N
    for j = 1:i
        K(i, j) = h((i-j)*dx)*dx;
    end
end





%% Local regularization of Volterra equations

delta_1 = sum(h(xgrid(1:2)))*dx/2;
prony_h_grid = h(xgrid)';
prony_g_grid = g(xgrid)';
A = tril(toeplitz(prony_h_grid(1:end-1)+prony_h_grid(2:end)))/2;


r = 10;
all_s       = zeros(r, 1);
all_A       = cell(r, 1);
all_f       = cell(r, 1);
delta       = zeros(r, 1);

for i = 1:r
    delta(i) = sum(prony_h_grid(1:i) + prony_h_grid(2:i+1))/2*dx;
    all_s(i) = delta(i) / delta(1);
    all_A{i} = tril(toeplitz(prony_h_grid(i:N-1-r+i)+prony_h_grid(i+1:N-r+i)))/2*dx;
    all_f{i} = -prony_g_grid(i:N-r+i-1);
end

AA = zeros(N-r, N-r);
ff = zeros(N-r, 1);
for i = 1:r
    if i == 1
        delta_i_1 = 0;
    else
        delta_i_1 = delta(i-1);
    end
    AA = AA + all_s(i)*(all_A{i} + delta_i_1*eye(N-r, N-r));
    ff = ff + all_s(i)*all_f{i};
end
    
% [~, cc] = L_curve(AA'*AA,AA'*ff,'RKHS',1);

cc = lsqminnorm(AA, ff);
%% Traditional regularization

g_grid = g(xgrid)';

A = K'*K;
b = -K'*g_grid;


c = A\b;
lambda = 1e-2;
c_RKHS = (A+lambda*eye(N))\b;


%%
figure;
subplot(121);hold on;
plot(xgrid, c, 'LineWidth', 3)
plot(xgrid, c_RKHS, 'LineWidth', 3)
plot(xgrid, result_gamma.gamma(xgrid), 'LineWidth', 3)
plot(xgrid, I.gamma(xgrid), 'LineWidth', 3)
plot(xgrid(1:N-r), cc, 'LineWidth', 3)
legend('LS','RKHS','Prony','true', 'reg Volterra')

subplot(122);hold on;
plot(log10(xgrid), log10(c), 'LineWidth', 3);
plot(log10(xgrid), log10(c_RKHS), 'LineWidth', 3)
plot(log10(xgrid), log10(result_gamma.gamma(xgrid)), 'LineWidth', 3)
plot(log10(xgrid), log10(I.gamma(xgrid)), 'LineWidth', 3)
plot(log10(xgrid(1:N-r)), log10(cc), 'LineWidth', 3)
legend('LS','RKHS','Prony','true','reg Volterra')

% subplot(122);hold on;
% loglog((c), 'LineWidth', 3);
% loglog((c_RKHS))
% loglog((result_gamma.gamma(xgrid)))
% loglog((I.gamma(xgrid)))
% loglog((cc))
% legend('LS','RKHS','Prony')


