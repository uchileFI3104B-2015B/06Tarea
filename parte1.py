from __future__ import division
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.sparse import (spdiags, coo_matrix, csc_matrix,
dia_matrix, dok_matrix, identity)


def solucion_inicial(x, t=0):
    return np.exp(-x**2 / 0.1)


def M(x):
    # me devuelve una matriz de len(x) x len(x), tridiagonal
    # con secuencia 1, -2, 1.
    ones = np.ones(len(x))
    M = spdiags([ones, -2*ones, ones], (-1, 0, 1), len(x), len(x)).tocsc()
    return M


def Imod(x):
    data = np.ones(len(x))
    data[0] = 0
    data[-1] = 0
    offset = 0
    shape = (len(x), len(x))
    Imod = dia_matrix((data, offset), shape=shape).todok()
    return Imod.tocsc()


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

# Matrices and vectores
M = M(x)                # Matriz de coeficientes
Im = Imod(x)            # matriz identidad modificada
I = identity(len(x))    # matriz identidad
