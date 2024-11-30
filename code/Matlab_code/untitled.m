E = trueInfo.E_true;
L = trueInfo.L_true;

E = (E + E')/2;
r = trueInfo.r_true;
global n
n = sysInfo.n;



[u, s, v] = svd(E);
ss = diag(s);
ss(1) = -ss(1);

U = u(:, 1:r)*diag(sqrt(ss(1:r)));




a = u*sqrt(s);
b = v*sqrt(s);
a = a(:, 1:r);
b = b(:, 1:r);

norm(a*b' - E)

%%
[uu, ss] = eig(E);
U = uu*sqrt(ss);
V = uu*conj(sqrt(ss));

U(:, 2:end-r+1) = [];
V(:, 2:end-r+1) = [];

norm(U*V' - E, 'fro')
norm(U'*U - V'*V, 'fro')


bk(U, 1)*bk(V, 2)'
BK(E, 1, 2)



%%
rho = all_pair_data.prony;
rho_0 = rho(:, :, 1, 1);
rho_1 = rho(:, :, 2, 1);



clc
% sum(conj(rho_0) .* BK(E, 1, 3)', 'all')
sum((rho_0) .* conj(BK(E, 1, 3)), 'all')
% sum((rho_0) .* (BK(E, 1, 3)), 'all')
% sum((rho_0)' .* conj(BK(E, 1, 3)), 'all')
trace(rho_0*BK(E, 1, 3)')


rho_1(1, 3)

% BK(E, 1, 1)
% rho_0
% 
% L(1, :)
% vec(rho_0)
% 
% 
% L(1, :) * vec(rho_0)
% 
% rho_1
%%
function Uj = bk(U, j)
    global n
    Uj = U((j-1)*n+1:j*n, :);
end

function Eij = BK(E, i, j)
    global n
    Eij = E((i-1)*n+1:i*n, (j-1)*n+1:j*n);
end