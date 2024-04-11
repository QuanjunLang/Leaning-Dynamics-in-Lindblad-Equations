function prony_result = Prony_fit_rho(all_rho_obs, obsInfo)
% Use prony method to fit the trajectory of each entry of rho

%% 
[n, ~, ~, M]    = size(all_rho_obs);
dt              = obsInfo.dt;
obs_len         = obsInfo.obs_len;
prony_p         = obs_len/2;
prony_p         = 5;
%% prony fitting of the density trajectories
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

end



