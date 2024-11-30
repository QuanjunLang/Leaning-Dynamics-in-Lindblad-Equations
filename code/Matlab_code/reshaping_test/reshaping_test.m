% % Example matrix A of size 4x4 (as a 2x2 block structure)
% K = rand(4, 4) + 1i*rand(4, 4);
% u = rand(4, 1) + 1i*rand(4, 1);
% v = rand(4, 1) + 1i*rand(4, 1);
% 
% K = u*v';

%%
A = rand(2,2) + 1i*rand(2,2);
A = A + A';
B = rand(2,2) + 1i*rand(2,2);
B = B + B';
% B = conj(A);
K = kron(A, B);


% Apply reshaping operator
clc
K
R_K = reshape_operator_col_first(K)
% vec_K = vec(A.')*vec(B.').';
vec_K = vec(A)*(vec(B).')
% vec_K = vec(B)*(vec(A).')
norm(R_K - vec_K)

RR_K = reshape_operator_col_first(R_K.')
K
norm(RR_K - K.')

% RR_K = reshape_operator(R_K);
% RR_K

%%
A = rand(2,2) + 1i*rand(2,2);
% A = A + A';
B = rand(2,2) + 1i*rand(2,2);
% B = B + B';

K = kron(A, B);


% Apply reshaping operator
clc
K
R_K = reshape_operator_row_first(K)
vec_K = vec(A.')*vec(B.').'
% vec_K = vec(A)*(vec(B).')

norm(R_K - vec_K)

RR_K = reshape_operator_row_first(R_K)
K
norm(RR_K - K)

% RR_K = reshape_operator(R_K);
% RR_K



%%
L = trueInfo.L_true;
E = trueInfo.E_true;

L
E
R_L = reshape_operator_col_first(L);
RR_L = reshape_operator_col_first(R_L);

norm(R_L - E)
norm(RR_L - L)


R_L_row = reshape_operator_row_first(L);
RR_L_row = reshape_operator_row_first(R_L_row);

norm(R_L_row - E)
norm(RR_L_row - L)