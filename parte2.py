#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve la ecuacion de Newell-Whitehead-Segel
La ecuación a resover es:
    dn/dt = gamma * d2n/dx2 + mu * n - mu * n ** 2
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# semilla
# np.random.seed(123)
np.random.seed(8)
# np.random.seed(888888888)
# funciones base
'''
def cond_ini(n, M):
    #establece las condiciones iniciales del problema, no se logra implementar
    for i in range(M):
        n[i] = np.random.uniform(low=-0.3, high=0.3, size=M)
        n[0] = 0 # cond iniciales
        n[-1] = 0
    pass
'''


def calc_b(b, M, r, h, mu):
    '''
    calcula el valor de b para el metodo de crank nicolson
    '''
    for j in range(1, M - 1):
        # se define a para que la linea siguiente no quede tan grande porque
        # no pude lograr que funcionara el corte de la linea :s
        a = h * mu * (n[j] - n[j] ** 3)
        b[j] = r * n[j+1] + (1 - 2 * r) * n[j] + r * n[j-1] + a
    pass


def calc_alpha_beta(alpha, beta, b, r, M):
    '''
    obtiene los valores de alpha y beta para el metodo de crank nicolson
    '''
    Amas = - 1 * r
    Acero = (1 + 2 * r)
    Amenos = - 1 * r
    alpha[0] = 0
    beta[0] = 0  # condiciones iniciales
    for i in range(1, M):
        alpha[i] = - Amenos / (Acero + Amenos * alpha[i-1])
        beta[i] = (b[i] - Amenos * beta[i-1]) / (Amenos * alpha[i-1] + Acero)
    pass


def avanza_tiempo(n, n_sig, alpha, beta, M):
    '''
    obtiene el comportamiento de la funcion al tiempo siguiente
    '''
    n_sig[0] = 0  # condiciones iniciales
    n_sig[-1] = 0
    for i in range(M - 12, 0, -1):
        n_sig[i] = alpha[i] * n_sig[i+1] + beta[i]
    pass
# main
# inicializacion
mu = 1.5
gamma = 0.001
dx = 0.002
dt = 0.002
N_x = 500
N_t = 3000
print N_x
print N_t
r = (gamma * dt) / (2 * dx ** 2)
n = np.random.uniform(low=-0.3, high=0.3, size=N_x)
n[0] = 0
n[-1] = 0
n_sig = np.zeros(N_x)
b = np.zeros(N_x)
alpha = np.zeros(N_x)
beta = np.zeros(N_x)
n_sol = np.zeros((N_t, N_x))
n_sol[0, :] = n.copy()
# iteracion
for i in range(1, N_t):
    calc_b(b, N_x, r, dt, mu)
    calc_alpha_beta(alpha, beta, b, r, N_x)
    avanza_tiempo(n, n_sig, alpha, beta, N_x)
    n = n_sig.copy()
    n_sol[i, :] = n.copy()

# plots
x = np.linspace(0, 1, N_x)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_t, 500):
    ax.plot(x, n_sol[i, :], "p")
    ax.plot(x, n_sol[i, :], "-p", label="t="+str(i*dt))
    ax.legend(loc='center left', bbox_to_anchor=(1., 0.5))
plt.subplots_adjust(left=None, bottom=None, right=0.8, top=None, wspace=None,
                    hspace=None)
plt.xlabel("x")
plt.ylabel("n(x)")
# plt.savefig("NWS1.png")  # semilla 123
plt.savefig("NWS2.png")  # semilla 8
# plt.savefig("NSW3.png")  # semilla 888888888
plt.show()
plt.draw()
