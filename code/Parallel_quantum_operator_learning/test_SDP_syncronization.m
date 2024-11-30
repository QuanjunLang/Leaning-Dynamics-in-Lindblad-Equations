clc
close all
clear all
rng(10)


% Load the required parameters
% E: List of pairs (i,j) in the set \mathcal{E}
% R: Cell array where R{i,j} is the corresponding R_ij matrix
% d: Dimension of each G_ii block
% n: Number of blocks (size of G is n*d x n*d)

% Example inputs (replace with your data)
E = [1, 2; 2, 3; 1, 3]; % Edges (i, j)
d = 2;                  % Dimension of diagonal blocks
n = 3;                  % Number of blocks
R = cell(n, n);         % Cell array for R_ij
for i = 1:n
    for j = 1:n
        R{i,j} = zeros(d); % Initialize R_ij as zeros (replace with actual values)
    end
end
R{1,2} = [1, 0; 0, 1]; % Example R_ij values
R{2,3} = [0, 1; -1, 0];
R{1,3} = [1, -1; 1, 1];

% CVX optimization
cvx_begin sdp
    variable G(n*d, n*d) symmetric
    minimize( sum(arrayfun(@(k) norm(G((E(k,1)-1)*d+1:E(k,1)*d, (E(k,2)-1)*d+1:E(k,2)*d) ...
        - R{E(k,1), E(k,2)}, 'fro'), 1:size(E,1))) )
    subject to
        for i = 1:n
            G((i-1)*d+1:i*d, (i-1)*d+1:i*d) == eye(d); % G_ii = I_d
        end
        G >= 0; % G is positive semidefinite
cvx_end

% Display results
disp('Optimized G:');
disp(G);
