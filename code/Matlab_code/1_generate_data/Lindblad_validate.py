# %%
import numpy as np
from qutip import Qobj, mesolve, sigmax, Options
import qutip
import warnings
import time
import pandas as pd
warnings.filterwarnings('ignore')
np.random.seed(0)

# np.random.seed(1)
def test():
    start_time = time.time()

    # load time steps
    try:
        times = np.array(tgrid)
    except NameError:
        times = np.linspace(0, 1, 10000)

    # load Hamiltonians
    try:
        H = Qobj(Hamiltonian)
        N = H.shape[0]
    except NameError:
        N = 6
        H = qutip.rand_herm(N, density=0.75)/5

    
    # load Jump operators
    try:
        temp_jump = np.array(Jump)
        if len(temp_jump.shape) == 2:
            temp_jump = np.expand_dims(temp_jump, 2)

        C = []
        for i in range(P):
            C.append(Qobj(temp_jump[:, :, i]))

    except NameError:
        P = 2
        C = []
        for i in range(P):
            # temp = 0.2*qutip.rand_herm(N, density=0.75, dims=None, pos_def=False)
            # temp = 0.2*qutip.rand_dm(N, density=0.75)

            temp = np.random.normal(0, 1, size=(N, N)) + 1j*np.random.normal(0, 1, size=(N, N))
            temp[0, 0] = temp[0, 0] - np.trace(temp)

            C.append(qutip.Qobj(0.05*temp))
    
    # load initial values
    try:
        all_rho_initial = Qobj(all_rho_0)
        Num_traj = len(all_rho_initial)
    except NameError:
        try:
            Num_traj = int(M)
        except NameError:
            Num_traj = 10

        all_rho_initial = []
        for m in range(Num_traj):
            all_rho_initial.append(qutip.rand_dm(N, density=0.8, seed=0))

    # generate trajectories
    all_traj = []
    for m in range(Num_traj):
        rho_initial = all_rho_initial[m]
        curr_traj = mesolve(H, rho_initial, times, C, options = Options(store_states=True))
        temp = []
        for item in curr_traj.states:
            temp.append(item.full())
        all_traj.append(temp)

    end_time = time.time()
    elapsed_time = end_time - start_time

    all_traj = np.array(all_traj)
    return all_traj, elapsed_time


# %%
a = test()

