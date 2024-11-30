ind = 3;


dt = sysInfo.dt;
n = sysInfo.n;

L = trueInfo.L_true;            % True Lindbladian
C = expm(trueInfo.L_true*dt);   % True Channel

% Using density matrix 
% rho0 = all_rho(:, :, 1, ind);
% rho1 = all_rho(:, :, 2, ind);
% rho0_prime = (rho1 - rho0)/dt;

% Using arbitrary matrix to test
rho0 = randn(n, n) + 1i*randn(n, n);
rho1 = reshape(C*vec(rho0), [n, n]);
rho0_prime = (rho1 - rho0)/dt;

%%
norm(C*vec(rho0) - vec(rho1))
norm(L*vec(rho0) - vec(rho0_prime))

%% Observable test

% Using position observable
% observable_ind_i = 1;
% observable_ind_j = 3;
% O = zeros(sysInfo.n, sysInfo.n);
% O(observable_ind_i, observable_ind_j) = 1;

% Using arbitrary observable
O = randn(n, n) + 1i*randn(n, n);



u = vec(O);
v = vec(rho0);


u' * C * v - u' * vec(rho1)
u' * L * v - u' * vec(rho0_prime)

%%
norm(L - L')    % L is not Hermitian
norm(C - C')    % C is not Hermitian

%% Trace Inner product
UV = u*v';
VU = v*u';          % UV = VU'


trace(C*VU) - u' * C * v        % u' * C * v = tr(C * (v*u'))
trace(L*VU) - u' * L * v

trace(UV'*C) - u' * C * v
trace(UV'*L) - u' * L * v

trace(VU*C) - u' * C * v
trace(VU*L) - u' * L * v

%%
trace(UV'*C) - vec(UV)'*vec(C)


%% Frobenius Inner product
sum(conj(UV).*C, 'all') - u' * C * v    % u'*C*v = <UV, C> = 
sum(conj(UV).*L, 'all') - u' * L * v

%% rearrangement
RUV = rearrangement_R(UV);
RC = rearrangement_R(C);
RL = rearrangement_R(L);

sum(conj(RUV).*RC, 'all') - u' * C * v
sum(conj(RUV).*RL, 'all') - u' * L * v

trace(RUV'*RC) - u' * C * v
trace(RUV'*RL) - u' * L * v

%%

norm(RUV - kron(O.', rho0'), 'fro')
norm(rearrangement_R(RUV) - UV, 'fro')
norm(rearrangement_R(rearrangement_R(RUV)) - rearrangement_R(UV), 'fro')


%%
norm(rearrangement_R(kron(O, rho0)) - vec(O.')*vec(rho0')', 'fro')
norm(rearrangement_R(kron(O, rho0)) - vec(O.')*vec(rho0.').', 'fro')
%%
sum(conj(kron(O.', rho0')).*RC, 'all') - u' * C * v

trace(kron(O.', rho0')'*RC) - u'*C*v