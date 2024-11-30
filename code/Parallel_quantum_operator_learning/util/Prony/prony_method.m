function result = prony_method(I)

x = I.obs_h_grid;
p = I.prony_p;
N = I.prony_N;
Ts = I.obs_dx;
% rho = I.rho;
drop_0 = I.drop_0;
polycoef_method = I.polycoef_method;
weight_method = I.weight_method;

root_normalization = I.root_normalization;
lambda_augmentation = I.lambda_augmentation;

if drop_0
    x_short = x(2:end);
    T = toeplitz(x_short(p:N-1),x_short(p:-1:1));
else
    T = toeplitz(x(p:N-1),x(p:-1:1));
end
switch polycoef_method
    case {'LS', 'RKHS', 'ID'}
        switch polycoef_method
            case {'LS'}
                a = -T\x(p+1:N);
            case {'RKHS', 'ID'}
                [~, a] = L_curve(T'*T, -T'*x(p+1:N), polycoef_method, 0);
        end
        
        % check for indeterminate forms
        indeterminate_form = sum(isnan(a) | isinf(a));
        if (indeterminate_form)
            disp('Coefficients of polynomial is nan')
            return;
        end
        
        % step 2
        c = transpose([1; a]);
        r = roots(c);
        
        
        
    case {'MP'}
        N = length(x);
        Y = hankel(x(1:end-p), x(end-p:end));
        
        Y1 = Y(:, 1:end-1);
        Y2 = Y(:, 2:end);
        
        r = eig(pinv(Y1)*Y2);
end

sigma = 1e-6;
% step 2.5
if root_normalization
    n = length(r);
    rr = r;
    for i = 1:n
        if abs(r(i)) >=1
            rr(i) = r(i)/(abs(r(i)) + sigma);
        end
    end
    r = rr;
end
plotON = 0;
% plot roots
if plotON

figure;
scatter(real(r), imag(r), 130, '.'); hold on;circle(0, 0, 1)
xrange = 1.4;
xlim([-xrange, xrange]);ylim([-xrange, xrange])
title('Prony modes')

end
%% step 3




% hold on;circle(0, 0, 1)
% xrange = 1.4;
% xlim([-xrange, xrange]);ylim([-xrange, xrange])
% title('Prony modes')

lam_straight = log(r)/Ts;

%
% figure;
% scatter(real(lam), imag(lam), 130, '.');

%%

% p = length(lam);


if drop_0
    Z = zeros(N-1, p);
    for i = 1:length(r)
        Z(:,i) = transpose(r(i) .^(1:N-1));
    end
    x = x_short;
else
    Z = zeros(N, p);
    for i = 1:length(r)
        Z(:,i) = transpose(r(i) .^(0:N-1));
    end
end

% Z = zeros(N, p);
% for i = 1:length(r)
%     Z(:,i) = transpose(exp(lam(i)*Ts) .^(0:N-1));
% end
% > lsalin(Z, x, [], [],

switch weight_method
    case {'LS'}
        % w = Z\x ;
        w = lsqminnorm(Z, x);
    case {'LS_freq'}
        D = diag((imag(lam_straight)).^2) + eye(p, p);
        A = Z'*Z;
        b = Z'*x;
        lambda = 1e-3;
        w = (A + lambda*D)\b;
        %         [~, h] = L_curve_standard_form(D\A, D\b, 1);
    case {'LS_rho'}
        rho_grid = rho(linspace(0, Ts, N)');
        ZZ = rho_grid.*Z;
        xx = rho_grid.*x;
        
        w = ZZ\xx;
        %         [~, w] = L_curve(ZZ'*ZZ, -ZZ'*xx, 'RKHS', 0);
    case {'LS_h0'}
        A = [Z;1000*lam_straight'];
        b = [x; 0];
        w = (A'*A)\(A'*b);
        %         h = lsqlin(Z,x(1:N),[],[],log(conj(r')),0, [], [], []);
        
    case {'LS_h0_new'}
        C = [-transpose(lam_straight(2:end))./lam_straight(1);eye(p-1)];
        A = Z*C;
        b = x;
        w_short = (A'*A)\(A'*b);
%         w_short = L_curve(A'*A, A'*b, 'ID', 1);
        %         h = lsqlin(Z,x(1:N),[],[],log(conj(r')),0, [], [], []);
        w = C*w_short;
        
        a = 1;
    case {'RKHS', 'ID'}
        [~, w] = L_curve(Z'*Z, Z'*x(1:N), weight_method, 0);
    case {'RKHS_h0', 'ID_h0'}
        A = [Z;100*log(conj(r'))];
        b = [x(1:N); 0];
        method = weight_method(1:end-3);
        [~, w] = L_curve(A'*A, A'*b, method, 1);
end




if lambda_augmentation
    n = length(r);
    r_bad = [];
    r_good = [];
    w_good = [];
    w_bad = [];
    
    for i = 1:n
        if real(r(i)) < 0 && imag(r(i)) == 0
            r_bad = [r_bad, r(i)];
            w_bad = [w_bad, w(i)];
        else
            r_good = [r_good, r(i)];
            w_good = [w_good, w(i)];
        end
    end
    lam = log(r_good)/Ts;
    new_w = w_good;
    
    
    for i = 1:length(r_bad)
        rr = r_bad(i);
        ww = w_bad(i);
        %         assert(abs(imag(ww)) < eps);
        lam = [lam, log(rr)/Ts];
        lam = [lam, conj(log(rr)/Ts)];
        new_w = [new_w, ww/2, ww/2];
    end
    
    %     old_w = w;
    %     old_lam = log(r)/Ts;
    
    w = transpose(new_w);
    lam = transpose(lam);
    p = length(w);
    
else
    lam = log(r)/Ts;
end


h = @(x) 0*x;
for t = 1:p
    h = @(x) h(x) + w(t).*exp(lam(t)*x);
end

% old_h = @(x) 0*x;
% for t = 3:length(old_w)
%     old_h = @(x) old_h(x) + old_w(t).*exp(old_lam(t)*x);
% end

h = @(x) h(x);


h_lap = @(x) 0*x;
for t = 1:p
    h_lap = @(x) h_lap(x) + w(t)./(x - lam(t));
end


% figure;hold on;
% xgrid = linspace(0, 10, 1000);
% plot(xgrid, h(xgrid));
% plot(xgrid, old_h(xgrid), 'o');
%
% [w(1:8), lam(1:8)]
% [old_w(3:end), old_lam(3:end)]


result.lam = lam;
result.r = r;
result.h = h;
result.w = w;
% result.a = a;
% result.c = c;
result.h_lap = h_lap;
result.h = h;
%
%
%
dh = @(x) 0*x;
for t = 1:p
    dh = @(x) dh(x) + lam(t).*w(t).*exp(lam(t)*x);
end
% dh = @(x) real(dh(x));

ddh = @(x) 0*x;
for t = 1:p
    ddh = @(x) ddh(x) + lam(t).*lam(t).*w(t).*exp(lam(t)*x);
end
% ddh = @(x) real(ddh(x));



result.dh = dh;
result.ddh = ddh;
%
%
plotON = 0;
%%
if plotON
    
    figure;
    subplot(121);hold on;
    plot(I.xgrid, I.h_grid, '-', 'LineWidth', 5, 'DisplayName', 'true')
    plot(I.obs_xgrid, I.obs_h_grid, '.', 'MarkerSize', 10, 'LineWidth', 3, 'DisplayName', 'obs')
    plot(I.xgrid, real(h(I.xgrid)), '-.', 'MarkerSize', 2, 'LineWidth', 2, 'DisplayName', 'est');
    legend;
    
    
    
    subplot(122);hold on;
    scatter(real(r), imag(r), 130, '.'); hold on;circle(0, 0, 1)
    xrange = 1.4;
    xlim([-xrange, xrange]);ylim([-xrange, xrange])
    title('poly')
    
end
end
