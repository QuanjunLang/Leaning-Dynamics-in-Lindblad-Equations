# %%
import numpy as np
from qutip import Qobj, mesolve, sigmax, Options
import qutip
import warnings
import time
import pandas as pd
warnings.filterwarnings('ignore')
rng = np.random.default_rng(12345)

# np.random.seed(1)
def test():
    start_time = time.time()

    try:
        N = int(n)
    except NameError:
        N = 4

    try:
        Num_jump = int(p)
    except NameError:
        Num_jump = 2
    
    try:
        Num_traj = int(M)
    except NameError:
        Num_traj = 3

    try:
        Num_observable = int(N_o)
    except NameError:
        Num_observable = 5

    try:
        times = np.array(tgrid)
    except NameError:
        times = np.linspace(0, 1, 1000)

    # print(N_o)


    # Hamiltonian
    H = qutip.rand_herm(N, density=0.75, seed=rng)/5
    
    # Jump operators
    C = []
    for i in range(Num_jump):
        temp = np.random.normal(0, 1, size=(N, N)) + 1j*np.random.normal(0, 1, size=(N, N))
        temp[0, 0] = temp[0, 0] - np.trace(temp) # make jump operators traceless
        C.append(qutip.Qobj(0.5*temp))

    # Observables
    O = []
    for m in range(Num_traj):
        Om = []
        for i in range(Num_observable):
            Om.append(qutip.rand_herm(N, seed=rng, density = 0.8))
        O.append(Om)


    end_time = time.time()
    elapsed_time = end_time - start_time

    all_traj = []
    all_initial = []
    for m in range(Num_traj):
        rho_initial = qutip.rand_dm(N, density=0.5, seed=rng) # initial states are randomly generated
        curr_traj = mesolve(H, rho_initial, times, C, e_ops = O[m])
        # curr_traj = mesolve(H, rho_initial, times, C)
        # temp = []
        # for item in curr_traj.states:
        #     temp.append(item.full())
        all_traj.append(np.array(curr_traj.expect))
        all_initial.append(rho_initial.full())


    # curr_traj.expect[0][0]
    # np.sum(np.transpose(O[-1][0].full())*rho_initial.full())

    # store trajectory
    all_traj = np.array(all_traj)

    # store Hamiltonian
    HH = H.full()
    HH = np.array(HH)

    # store Jumps
    CC = []
    for item in C:
        CC.append(item.full())
    # CC = np.array(CC)

    # store Observables
    OO = []
    for Om in O:
        OOm = []
        for item in Om:
            OOm.append(item.full())
        OO.append(OOm)
    # OO = np.array(OO)
    
    return all_traj, HH, CC, OO, all_initial, elapsed_time


a = test()


# %%
