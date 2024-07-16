function sysInfo = update_sys(sysInfo)


sysInfo.L       = sysInfo.steps+1;
sysInfo.T       = sysInfo.steps*sysInfo.dt;
sysInfo.tgrid   = 0:sysInfo.dt:sysInfo.steps*sysInfo.dt;

sysInfo.FULL_STATE = 1;
if isfield(sysInfo, 'N_o')
    if sysInfo.N_o < sysInfo.n^2
        sysInfo.FULL_STATE = 0;
    end
else
    sysInfo.N_o = sysInfo.n^2;
end


end