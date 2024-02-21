# %%
import matplotlib.pyplot as plt
import numpy as np
from qutip import Bloch, about, basis, mesolve, sesolve, sigmam, sigmax, sigmay, sigmaz

import matplotlib.pyplot as plt
# import matplotlib as mpl
# from matplotlib import cm

from qutip import *
from qutip.piqs import *


# %%Example of Schodinger equation
H = 2*np.pi * 0.1 * sigmax()

h = np.array([[0, 1], [1, 0]])
H = Qobj(h)

times = np.linspace(0.0, 1.0, 2)
rho0 = Qobj([[1,3],[3, 2]])


result_0 = mesolve(H, rho0, times, options = Options(store_states=True))



# %%
