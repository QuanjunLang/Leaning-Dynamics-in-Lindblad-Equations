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
        N = 20

    try:
        Num_jump = int(p)
    except NameError:
        Num_jump = 2
    
    try:
        Num_traj = int(M)
    except NameError:
        Num_traj = 200

    try:
        times = np.array(tgrid)
    except NameError:
        times = np.linspace(0, 0.00001, 10)

    # try:
    #     Diff_Obs_Flag = UseDiffObs
    # except NameError:
    #     Diff_Obs_Flag = 0


    # Hamiltonian
    H = qutip.rand_herm(N, density=1, seed=rng)/5
    
    # Jump operators
    C = []
    for i in range(Num_jump):
        temp = np.random.normal(0, 1, size=(N, N)) + 1j*np.random.normal(0, 1, size=(N, N))
        temp[0, 0] = temp[0, 0] - np.trace(temp) # make jump operators traceless
        # temp = qutip.rand_herm(N, density=1, seed=rng)/5
        C.append(qutip.Qobj(0.2*temp))

    
# Using different observables to observe a given trajectory multiple times
    print('Generate data: first row, col, diagonal observables for all trajs')
    # Observables
    O_R = []
    O_I = []
    O_D = []
    for i in range(N):
        O_D.append(generate_observable_diag(N, i))
    for i in range(1, N):
        O_R.append(generate_observable_sym(N, 0, i))
        O_I.append(generate_observable_asym(N, 0, i))
    O = O_D + O_R + O_I

    all_traj = []
    all_initial = []
    for m in range(Num_traj):
        rho_initial = qutip.rand_dm(N, density=0.9, seed=rng) # initial states are randomly generated
        curr_traj = mesolve(H, rho_initial, times, C, e_ops = O)
        all_traj.append(np.array(curr_traj.expect))
        all_initial.append(rho_initial.full())

    # store Observables
    OO = []
    for item in O:
        OO.append(item.full())

    

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



    end_time = time.time()
    elapsed_time = end_time - start_time



    return all_traj, HH, CC, all_initial, elapsed_time, OO


def generate_observable_sym(N, i, j):
    """
    Generate the quantum observable E_ij + E_ji for a Hilbert space of dimension N.

    Parameters:
        N (int): Dimension of the Hilbert space.
        i (int): Row index (0-based).
        j (int): Column index (0-based).

    Returns:
        Qobj: The observable E_ij + E_ji.
    """
    # Basis vectors
    ket_i = qutip.basis(N, i)
    ket_j = qutip.basis(N, j)
    
    # E_ij = |i><j| and E_ji = |j><i|
    E_ij = ket_i * ket_j.dag()
    E_ji = ket_j * ket_i.dag()
    
    # Symmetric sum
    observable = E_ij + E_ji
    return observable

def generate_observable_asym(N, k, l):
    """
    Generate the quantum observable iE_kl - iE_lk for a Hilbert space of dimension N.

    Parameters:
        N (int): Dimension of the Hilbert space.
        k (int): Row index (0-based).
        l (int): Column index (0-based).

    Returns:
        Qobj: The observable iE_kl - iE_lk.
    """
    # Basis vectors
    ket_k = qutip.basis(N, k)
    ket_l = qutip.basis(N, l)
    
    # E_kl = |k><l| and E_lk = |l><k|
    E_kl = ket_k * ket_l.dag()
    E_lk = ket_l * ket_k.dag()
    
    # Antisymmetric observable
    observable = 1j * E_kl - 1j * E_lk
    return observable

def generate_observable_diag(N, k):
    """
    Generate the quantum observable E_kk for a Hilbert space of dimension N.

    Parameters:
        N (int): Dimension of the Hilbert space.
        k (int): Index for the diagonal element (0-based).

    Returns:
        Qobj: The observable E_kk.
    """
    # Basis vector |k>
    ket_k = qutip.basis(N, k)
    
    # Projector E_kk = |k><k|
    E_kk = ket_k * ket_k.dag()
    return E_kk


a = test()


# %%
