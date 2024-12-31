# %%
import numpy as np
import qutip
import warnings
import time
warnings.filterwarnings('ignore')
# rng = np.random.default_rng(12345)

# np.random.seed(1)
def test():
    start_time = time.time()

    try:
        M = int(m)
    except NameError:
        M = 17

    try:
        N = int(n)
    except NameError:
        N = 8

    try:
        seed = int(random_seed)
    except NameError:
        seed = 1


    O = [qutip.rand_herm(N, seed=seed+_).full() for _ in range(M)]
    rho0 = [qutip.rand_dm(N, density=1, seed=seed+_).full() for _ in range(M)]

    # Observables
  
    return [np.array(O), np.array(rho0)]





# %%
a = test()