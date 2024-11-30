function R_A = reshape_operator_col_first(A)
    % A is an (N^2 x N^2) matrix where N = n * n
    % A is divided into n x n blocks, and each block A_ij is vectorized
    %
    % Input: A - a matrix of size (N^2 x N^2) where N = n^2
    %        n - blocking parameter, which gives us the number of blocks
    % Output: R_A - rearranged matrix where each block A_ij is vectorized

    % Get the size of A
    [N2, ~] = size(A);
    
    % N is the square root of N2, assuming A is (N^2 x N^2)
    N = sqrt(N2);
    
    % Preallocate the reshaped matrix R_A
    R_A = zeros(N2, N^2);  % N^2 x n matrix to hold vec(A_ij)
    
    % Loop over each block of size n x n
    for j = 1:N
        % Collect the j-th block column
        for i = 1:N
            A_ij = A((i-1)*N+1:i*N, (j-1)*N+1:j*N);  % Extract block A_ij
            R_A((j-1)*N+i, :) = reshape(A_ij, 1, []);
        end
        
    end
end
