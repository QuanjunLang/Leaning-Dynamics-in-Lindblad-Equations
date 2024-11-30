clear all
clc
close all

%%

N = 4;
Nk = 2;
Nm = 1;
Np = Nk - Nm;


%% Prepare some basis vectors
A = rand(N^2, N^2) + 1i*rand(N^2, N^2);
[W_temp, ~, ~] = svd(A);
W = W_temp(:, 1:Nk);


if Nm == 0
    E = diag(ones(Np, 1));
else
    E = diag([-1*ones(Nm, 1); ones(Np, 1) ]);
end

K = W*E*W';
K = (K+K')/2;

%% blocks 
K_block = cell(N, N);

for i = 1:N
    for j = 1:N
        K_block{i, j} = K((i-1)*N+1:i*N, (j-1)*N+1:j*N);
    end
end




%% 
W_block = cell(N, 1);
for i = 1:N
    W_block{i} = W((i-1)*N+1:i*N, :);
end


%% check equivalence 
for i = 1:N
    for j = 1:N
        K_ij = W_block{i} * E * W_block{j}';
        error(i, j) = norm(K_ij - K_block{i, j});
    end
end

error


%% From block to whole
WW = zeros(size(W));
WW_block = cell(N, 1);
for i = 1:N
    [U, S] = eig(K_block{i, i});
    s = diag(S);
    ind_m = 1:Nm;
    ind_p = N-Np+1:N;
    

    WW_block{i} = [U(:, ind_m)*sqrt(-1*diag(s(ind_m))), U(:, ind_p)*sqrt(diag(s(ind_p)))];
    WW((i-1)*N+1 : i*N, :) = WW_block{i};
end


KK = WW*E*WW';







% C_block = cell(N, 1);
% for i = 1:N
%     [C_block{i} = svd(W_block{i} * E * W_block{i}');
% end


%%
for i =1:N
    for j = 1:N
        a(i, j) = rank(K_block{i, j});
    end
end



%% Deterministic recomstruction
K_col_1 = K(:, 1:N);
[AA, S, BB] = svd(K_col_1);

KK_determine = zeros(N^2, N^2);
KK_determine(:, 1:N) = K(:, 1:N);

col_basis = AA(:, 1:Nk);
col_basis_block = cell(N, 1);
for i = 1:N
    col_basis_block{i} = AA((i-1)*N+1 : i*N, 1:Nk);
    rr(i) = rank(col_basis_block{i});
end



for i = 2:N
    coef{i} = col_basis_block{i} \ K_block{i, i};
    KK_determine(:, (i-1)*N+1:i*N) = col_basis * coef{i};
end

norm(KK_determine - K)
