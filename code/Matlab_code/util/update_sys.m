function sysInfo = update_sys(sysInfo)


sysInfo.L       = sysInfo.steps+1;
sysInfo.T       = sysInfo.steps*sysInfo.dt;
sysInfo.tgrid   = 0:sysInfo.dt:sysInfo.steps*sysInfo.dt;


if isfield('obs_std', sysInfo)
    obs_std = sysInfo.obs_std;
else
    obs_std = 0;
end


end