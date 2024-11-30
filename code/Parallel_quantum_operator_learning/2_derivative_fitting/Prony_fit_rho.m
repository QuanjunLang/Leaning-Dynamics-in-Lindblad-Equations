function prony_result = Prony_fit_rho(all_rho_obs, obsInfo)
% Use prony method to fit the trajectory of each entry of rho

if length(size(all_rho_obs)) == 2
    FULL_STATE = 0;
    UseDiffObs = 1;
    [~, M]    = size(all_rho_obs);
elseif length(size(all_rho_obs)) == 3
    FULL_STATE = 0;
    UseDiffObs = 0;
    [N_o, ~, M]    = size(all_rho_obs);
else
    FULL_STATE = 1;
    [n, ~, ~, M]    = size(all_rho_obs);
end

%% extract parameters

dt              = obsInfo.dt;
obs_len         = obsInfo.obs_len;
prony_p         = obs_len/2;
% prony_p         = 5;



%% prony fitting of the density trajectories
if FULL_STATE
    prony_result = cell(n, n, M);
    for m = 1:M
        rho = all_rho_obs(:, :, :, m);
        for ind_1 = 1:n
            for ind_2 = 1:n
                Prony_I.prony_p = prony_p;
                Prony_I.prony_N = obs_len;
                Prony_I.obs_dx  = dt;
                Prony_I.obs_h_grid = squeeze(rho(ind_1, ind_2, :));
                Prony_I.polycoef_method = 'MP';
                Prony_I.weight_method = 'LS';
                Prony_I.root_normalization = 0;
                Prony_I.lambda_augmentation = 0;
                Prony_I.drop_0 = 0;
                prony_result{ind_1, ind_2, m} = prony_method(Prony_I);
            end
        end
    end

elseif ~UseDiffObs
    prony_result = cell(N_o, M);
    for m = 1:M
        rho = all_rho_obs(:, :, m);
        for ind_o = 1:N_o
            Prony_I.prony_p = prony_p;
            Prony_I.prony_N = obs_len;
            Prony_I.obs_dx  = dt;
            Prony_I.obs_h_grid = squeeze(rho(ind_o, :))';
            Prony_I.polycoef_method = 'MP';
            Prony_I.weight_method = 'LS';
            Prony_I.root_normalization = 0;
            Prony_I.lambda_augmentation = 0;
            Prony_I.drop_0 = 0;

            prony_result{ind_o, m} = prony_method(Prony_I);
        end

    end

else
    prony_result = cell(M);
    for m = 1:M
        rho = all_rho_obs(:, m);
        Prony_I.prony_p = prony_p;
        Prony_I.prony_N = obs_len;
        Prony_I.obs_dx  = dt;
        Prony_I.obs_h_grid = rho;
        Prony_I.polycoef_method = 'MP';
        Prony_I.weight_method = 'LS';
        Prony_I.root_normalization = 0;
        Prony_I.lambda_augmentation = 0;
        Prony_I.drop_0 = 0;

        prony_result{m} = prony_method(Prony_I);
    end

end



end





