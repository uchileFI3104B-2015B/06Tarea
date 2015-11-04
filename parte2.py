#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve un problema simple de diffusion en 1D.
La ecuación a resover es:
    dn/dt = gama * d2n/dx2 + un - un^3;
    n(0,x) = np.random.uniform(low=-0.3, high=0.3, size=N_steps);
    n(t, 0) = 0 n(t, 1) = 0
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(386)


def inicializa_n(n, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        n[i] = semilla[i]
    n[0] = 0  # c de borde
    n[-1] = 0  # cndiciones de borde


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = (r * n[j+1] + (1-2*r) * n[j] + r * n[j-1] +
                dt * mu * (n[j] - n[j]**3))  # solucion


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alpha, beta, N_steps):
    n_next[0] = 0
    n_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]

# Main

# setup
gamma = 0.001
mu = 1.5
t_inicial = 0
t_final = 5
x_inicial = 0
x_final = 1

N_steps = 500
N_pasos_temporales = 1000

semilla = np.random.uniform(low=-0.3, high=0.3, size=N_steps)


dt = (t_final - t_inicial) / (N_pasos_temporales - 1)
h = (x_final - x_inicial) / (N_steps-1)
r = dt / (2 * h**2)
r = r * gamma
n = np.zeros(N_steps)
n_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_n(n, N_steps, h)

# Queremos guardar las soluciones en cada paso
n_solucion = np.zeros((N_pasos_temporales, N_steps))
n_solucion[0, :] = n.copy()

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(n, n_next, alpha, beta, N_steps)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()


# Plots


x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 10):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(-1, 1)
ax.set_xlabel('$\ x $', fontsize=15)
ax.set_ylabel('$\ n $', fontsize=15)
ax.set_title('$\ Curvas \ de \ n(t,x)$', fontsize=13)
plt.grid(True)
plt.show()
plt.draw()
