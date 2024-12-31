# %%
import numpy as np
from qutip import Qobj, mesolve, sigmax, Options, rand_dm, rand_super_bcsz, to_kraus, operator_to_vector, vector_to_operator, expect
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
        N = 8

    try:
        Kraus_rank = int(p)
    except NameError:
        Kraus_rank = 2
    
    try:
        Num_traj = int(M)
    except NameError:
        Num_traj = 300

    try:
        seed = int(random_seed)
    except NameError:
        seed = 1

    rng = np.random.default_rng(seed)

    # Super operator
    E = rand_super_bcsz(N, rank=Kraus_rank, seed=rng)
    K = qutip.to_kraus(E)
    Choi = qutip.to_choi(E)
    
    # Using different observables to observe a given trajectory multiple times
    print('Generate data: first row, col, diagonal observables for all trajs')
    # Observables
    O_D = [generate_observable_diag(N, 0)]
    O_R = [generate_observable_sym(N, 0, i) for i in range(1, N)]
    O_I = [generate_observable_asym(N, 0, i) for i in range(1, N)]
    O = O_D + O_R + O_I 

    # 
    rho_initial = [rand_dm(N, density=1, seed=rng) for _ in range(Num_traj)]
    rho_outputs = [vector_to_operator(E * operator_to_vector(rho_0)) for rho_0 in rho_initial]
    result = np.array([[expect(op, rho_1) for op in O] for rho_1 in rho_outputs])

    # t0 = time.time()
    # rho_0 = rho_initial[0]
    # rho_1 =vector_to_operator(E * operator_to_vector(rho_0))
    # expect(O[0], rho_1)
    # t1 = time.time()
    # t = t1 - t0
    # print(t)

    end_time = time.time()
    elapsed_time = end_time - start_time

    # return all_traj, H.full(), [c.full() for c in C], all_initial, elapsed_time, [o.full() for o in O]
    return result, E.full(), Choi.full(), [rho_0.full() for rho_0 in rho_initial], elapsed_time, [o.full() for o in O], np.array([rho1.full() for rho1 in rho_outputs]), [k.full() for k in K]


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
