clc
close all
clear all
rng(1)


%%
r = 2;
n = 2;
N = n^2;
all_A = zeros(n, n, r);
all_B = zeros(n, n, r);

for i = 1:r
    all_A(:, :, i) = rand(n, n);
    all_B(:, :, i) = rand(n, n);
end

Test1       = zeros(n^2, n^2);
Test1_all   = zeros(n^2, n^2, r);
for i = 1:r
    A = all_A(:, :, i);
    B = all_B(:, :, i);
    Test1_all(:, :, i) = vec(A)*vec(B)';
    Test1 = Test1 + Test1_all(:, :, i);
end

Test2_kron = zeros(n^2, n^2);
for i = 1:r
    A = all_A(:, :, i);
    B = all_B(:, :, i);
    Test2_kron = Test2_kron + kron(A, B);
end

Test2 = mat_vec(Test2_kron);



Test1 - Test2

%%
clc
close all
clear all

a = sym('a', [1,4]);
c = sym('c', [1,4]);
b = sym('b', [1,4]);

%%
A = [a(1), a(3); a(2), a(4)];
B = [b(1), b(3); b(2), b(4)];
C = [c(1), c(3); c(2), c(4)];

D = A*C*B



%%
R_mat = zeros(N, N);
for i = 1:N
    for j = 1:N
        temp = zeros(N, N);
        temp(i, j) = 1;
        mat_vec(temp)
    end
end






%%

function R = mat_vec(K)

N = length(K);
n = sqrt(N);

K_blocks = permute(reshape(K, n, n, n, n), [1, 3, 2, 4]);
R = permute(reshape(K_blocks, N, N), [2, 1]);

% R = zeros(N, N);
% k = 1;
% for i = 1:n
%     for j = 1:n
%         [i, j];
%         K_cur = K(((j-1)*n+1):j*n, ((i-1)*n+1):i*n);
%         R_cur = vec(K_cur);
%         R(k, :) = R_cur';
%         k = k+ 1;
%     end
% end
end





