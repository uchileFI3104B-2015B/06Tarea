from __future__ import division
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.sparse import (spdiags, coo_matrix, csc_matrix,
dia_matrix, dok_matrix, identity)



# Coeficientes y constantes
x = np.linspace(0, 1, 500)
x_range = (np.max(x) - np.min(x))
dx = x_range / len(x)
gamma = 0.001
mu = 1.5
s = gamma / dx**2
dt = 0.01
valor_izq = 1.0
valor_der = 0.0
