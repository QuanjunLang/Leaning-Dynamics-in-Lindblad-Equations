N = length(result_prony.E_est);
niter = length(result_prony.all_A);
all_UUT = zeros(N, N, niter);
all_VVT = zeros(N, N, niter);
all_E = zeros(N, N, niter);
   
r = trueInfo.r_true;
for i = 1:niter
    A_0 = result_prony.all_A(:, :, :, i);
    B_0 = result_prony.all_B(:, :, :, i);
    

    UU0 = zeros(N, N);
    VV0 = zeros(N, N);
    EE0 = zeros(N, N);

    for k = 1:r
        UU0 = UU0 + vec(A_0(:, :, k))*vec(permute(A_0(:, :, k), [2, 1])).';
        VV0 = VV0 + vec(B_0(:, :, k))*vec(permute(B_0(:, :, k), [2, 1])).';
        EE0 = EE0 + vec(B_0(:, :, k))*vec(permute(A_0(:, :, k), [2, 1])).';
    end

    all_UUT(:, :, i) = UU0;
    all_VVT(:, :, i) = VV0;
    all_E(:, :, i) = EE0;

    all_diff(i) = norm(UU0 - VV0, 'fro');

    
end
