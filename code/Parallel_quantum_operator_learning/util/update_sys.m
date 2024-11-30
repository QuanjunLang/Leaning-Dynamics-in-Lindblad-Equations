function sysInfo = update_sys(sysInfo)

sysInfo.L       = sysInfo.steps+1;
sysInfo.T       = sysInfo.steps*sysInfo.dt;
sysInfo.tgrid   = 0:sysInfo.dt:sysInfo.steps*sysInfo.dt;
sysInfo.channel_dt      = sysInfo.channel_dt_rate * sysInfo.dt;

switch sysInfo.observable_option
    case 'Full_state'
        sysInfo.N_o = sysInfo.n^2;
    case 'Multiple_random_observables'
        a = 1;
    case 'Single_random_observable'
        sysInfo.N_o = 1;
    case 'First_row_col_diag'
        sysInfo.N_o = 3*sysInfo.n - 2;
end



end