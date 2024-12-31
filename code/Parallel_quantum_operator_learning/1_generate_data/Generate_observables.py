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
        Num_O = int(N_o)
    except NameError:
        Num_O = 17

    try:
        N = int(n)
    except NameError:
        N = 8

    try:
        seed = int(random_seed)
    except NameError:
        seed = 1


    O = [qutip.rand_herm(N, seed=seed+_) for _ in range(Num_O)]

    # Observables
  
    return np.array([o.full() for o in O])

    return 1




# %%
a = test()