# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 23:29:51 2015

@author: splatt
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

'''
Este script resuelve un problema de modelamiento sobre el comportamiento de
una especie animal.
La ecuación a resolver es:
    dn/dt = Y * d2n/dx2 + u * n - u * (n²);
Con condiciones de borde:
n(t,0) = 1; n(t, 1) = 0; n(0,x) = e^(-x²/0.1);
'''


def inicializa_n(n, N_steps, h):
    '''
    Rellena n con las condiciones iniciales del problema.
    '''
    for i in range(N_steps):
        x = i * h
        n[i] = np.exp(-((x**2)/0.1))
    n[0] = 1
    n[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = (r * n[j+1] + (1-2*r) * n[j] +
                r * n[j-1] + dt*u*(n[j] - n[j]**2))


def calcula_alpha_y_beta(alpha, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde n(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alpha, beta, N_steps):
    n_next[0] = 1
    n_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]

# Main

# setup
N_steps = 501
h = 1 / (N_steps - 1)
dt = 0.01
N_pasos_temporales = 400
Y = 0.001  # gamma
u = 1.5  # mu
r = Y * dt / 2 / (h**2)

n = np.zeros(N_steps)
n_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_n(n, N_steps, h)


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

for i in range(0, N_pasos_temporales, 50):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(0, 1)
plt.xlabel('x (posicion)')
plt.ylabel('n(x,t) (densidad)')

fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax2.pcolormesh(X, Y, n_solucion)
plt.xlabel('x (posicion)')
plt.ylabel('t (tiempo)')

plt.show()
plt.draw()
