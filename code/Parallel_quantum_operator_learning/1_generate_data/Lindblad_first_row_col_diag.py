# %%
import numpy as np
from qutip import Qobj, mesolve, sigmax, Options, rand_dm
import qutip
import warnings
import time
from scipy.sparse import csr_matrix
import pandas as pd
from joblib import Parallel, delayed
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
        Num_traj = 300

    try:
        times = np.array(tgrid)
    except NameError:
        times = np.linspace(0, 0.00001, 2)

    # Hamiltonian
    H = qutip.rand_herm(N, density=0.5, seed=rng)/5
    
    # Jump operators
    C = []
    for i in range(Num_jump):
        temp = np.random.normal(0, 2, size=(N, N)) + 1j*np.random.normal(0, 2, size=(N, N))
        temp[0, 0] = temp[0, 0] - np.trace(temp) # make jump operators traceless
        # temp = qutip.rand_herm(N, density=1, seed=rng)/5
        C.append(qutip.Qobj(0.5*temp))

    
# Using different observables to observe a given trajectory multiple times
    print('Generate data: first row, col, diagonal observables for all trajs')
    # Observables
    O_D = [generate_observable_diag(N, i) for i in range(N)]
    O_R = [generate_observable_sym(N, 0, i) for i in range(1, N)]
    O_I = [generate_observable_asym(N, 0, i) for i in range(1, N)]
    O = O_D + O_R + O_I


    # Trajectory simulation
    all_traj = np.zeros((Num_traj, len(O), len(times)))
    all_initial = np.zeros((Num_traj, N, N), dtype=complex)
    options = Options(nsteps=1000, num_cpus=11)  # Parallel processing


    ####################################################################
    # Function to compute a single trajectory
    def compute_trajectory(m, H, N, times, C, O, rng_seed):
        # Generate initial state with the provided random seed
        rng = np.random.default_rng(rng_seed)
        rho_initial = rand_dm(N, density=0.9, seed=rng)  # Initial state
        
        # Solve the master equation for this trajectory
        curr_traj = mesolve(H, rho_initial, times, C, e_ops=O)
        
        # Return the trajectory and initial state
        return np.array(curr_traj.expect), rho_initial.full()

    # Number of trajectories
    # Num_traj = 300

    # Parallel computation of all trajectories
    all_results = Parallel(n_jobs=11)(  # Adjust n_jobs to the number of CPUs
        delayed(compute_trajectory)(m, H, N, times, C, O, 12345 + m)  # Unique seed for each trajectory
        for m in range(Num_traj)
    )

    # Unpack the results
    all_traj = np.array([result[0] for result in all_results])  # Expectation values
    all_initial = np.array([result[1] for result in all_results])  # Initial states
    ####################################################################


    # for m in range(Num_traj):
    #     rho_initial = rand_dm(N, density=0.9, seed=rng)  # Initial state
    #     # curr_traj = mesolve(H, rho_initial, times, C, e_ops=O, options=options)
    #     curr_traj = mesolve(H, rho_initial, times, C, e_ops=O)
    #     all_traj[m] = np.array(curr_traj.expect)
    #     all_initial[m] = rho_initial.full()

    end_time = time.time()
    elapsed_time = end_time - start_time

    return all_traj, H.full(), [c.full() for c in C], all_initial, elapsed_time, [o.full() for o in O]


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
    # # Convert to sparse explicitly (CSR format)
    # sparse_observable = Qobj(csr_matrix(observable.full()), dims=observable.dims, type='csr')
    # return sparse_observable

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
    # # Convert to sparse explicitly (CSR format)
    # sparse_observable = Qobj(csr_matrix(observable.full()), dims=observable.dims, type='csr')
    # return sparse_observable



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
    # Convert to sparse explicitly (CSR format)
    # sparse_observable = Qobj(csr_matrix(E_kk.full()), dims=E_kk.dims, type='csr')
    # return sparse_observable

    return E_kk


a = test()


# %%
