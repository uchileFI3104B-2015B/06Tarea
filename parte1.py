from __future__ import division
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.sparse import (spdiags, coo_matrix, csc_matrix,
dia_matrix, dok_matrix, identity)


def solucion_inicial(x, t=0):
    return np.exp(-x**2 / 0.1)


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

# Estado inicial
n_inicial = np.matrix(solucion_inicial(x)).reshape((len(x), 1))
n_inicial[0, 0] = valor_izq
n_inicial[-1, -1] = valor_der
n = n_inicial
