function plot_prony_rho_one_sample(sysInfo, obsInfo, all_rho_prony, all_rho, all_rho_obs)

FULL_STATE = sysInfo.FULL_STATE;

if FULL_STATE
    ind_1 = randi(sysInfo.n);
    ind_2 = randi(sysInfo.n);
    m = randi(sysInfo.M);

    rho_prony   = all_rho_prony{ind_1, ind_2, m};
    rho         = squeeze(all_rho(ind_1, ind_2, :, m));
    rho_obs     = squeeze(all_rho_obs(ind_1, ind_2, :, m));

else
    ind_o = randi(sysInfo.N_o);
    m = randi(sysInfo.M);

    rho_prony   = all_rho_prony{ind_o, m};
    rho         = squeeze(all_rho(ind_o, :, m));
    rho_obs     = squeeze(all_rho_obs(ind_o, :, m));
end

plot_prony_rho_new(rho_prony, rho, rho_obs, sysInfo, obsInfo);

end


