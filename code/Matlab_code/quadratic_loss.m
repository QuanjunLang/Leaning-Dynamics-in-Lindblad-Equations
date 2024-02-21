H   = sym('h', [2, 2]);
rho = sym('rho', [2, 2]);
A   = sym('a', [2, 2]);

%%
E = norm(A  - H*rho + rho*H,  'fro')^2;


%%
pretty(diff(E, H(1, 1)))