function prony = update_prony(prony)


obs_dx =    prony.obs_dx;
obs_std =   prony.obs_std;
p =         prony.p;
N =         prony.N;

prony.obs_xgrid = (0:obs_dx:(N-1)*obs_dx)';

end

