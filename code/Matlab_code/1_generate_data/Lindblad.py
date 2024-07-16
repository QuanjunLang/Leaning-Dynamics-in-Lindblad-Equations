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
    # Some defaul settings

    # s     = 0.01
    # H     = 2*np.pi * s * sigmax()
    # times = np.linspace(0.0, 1.0, 10)
    # rho0  = Qobj([[1,3],[3, 2]])
    # C     = [np.sqrt(0.05) * sigmax(), 0.1*sigmaz()]
    # H = Qobj(np.array(h))
    # print(N)
    # print(len(N))
    # rho_initial = Qobj(np.array(rho0))

    # CC = np.array(c)
    # C = []
    # for item in CC:
    #     C.append(Qobj(item))

    try:
        N = int(n)
    except NameError:
        N = 4

    try:
        P = int(p)
    except NameError:
        P = 2
    
    try:
        Num_traj = int(M)
    except NameError:
        Num_traj = 10


    try:
        times = np.array(tgrid)
    except NameError:
        times = np.linspace(0, 1, 10000)

    # N = int(n)
    # P = int(p)
    # Num_traj = int(M)
    # times = np.array(tgrid)

    # Hermitian term
    # H = qutip.rand_herm(N, density=0.75, dims=None, pos_def=False)/5
    H = qutip.rand_herm(N, density=0.75)/5
    
    # Lindbladian term
    C = []
    for i in range(P):
        # C.append(0.2*qutip.rand_herm(N, density=0.75, dims=None, pos_def=False))
        # C.append(0.2*qutip.rand_dm(N, density=0.75))

        temp = np.random.normal(0, 1, size=(N, N)) + 1j*np.random.normal(0, 1, size=(N, N))
        #
        temp[0, 0] = temp[0, 0] - np.trace(temp) 
        C.append(qutip.Qobj(0.05*temp))

        # C.append(0.2*qutip.destroy(N))
    
    # Generate Multiple trajectories

    # Initial state
    # all_traj = []
    # for m in range(Num_traj):
    #     rho_initial = qutip.rand_dm(N, density=0.5)
    #     curr_traj = mesolve(H, rho_initial, times, C, options = Options(store_states=True))
    #     all_traj.append(curr_traj.states)

    end_time = time.time()
    elapsed_time = end_time - start_time

    all_traj = []
    for m in range(Num_traj):
        rho_initial = qutip.rand_dm(N, density=0.5)
        curr_traj = mesolve(H, rho_initial, times, C, options = Options(store_states=True))
        temp = []
        for item in curr_traj.states:
            temp.append(item.full())
        all_traj.append(temp)

    all_traj = np.array(all_traj)
    H = H.full()
    CC = []
    for item in C:
        CC.append(item.full())
    # np.savetxt('data.csv', all_traj,  delimiter=',')
    return all_traj, H, CC, elapsed_time


# %%
a = test()


# %%
