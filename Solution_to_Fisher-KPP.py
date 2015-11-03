#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Este script es para resolver la ecuación de Fisher-KPP
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

#-implementar Crank-N

def inicializa_N(N, N_steps, h):
    '''
    Rellena N con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean uno y cero.
    '''
    for i in range(N_steps):
        x = i * h
        N[i] = np.exp((-x**2) / 0.1)
    N[0] = 1
    N[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * N[j+1] + (1-2*r) * N[j] + r * N[j-1]


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde T(t, 0) = 1 ¿?
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(N, N_next, alpha, beta, N_steps, mu=1.5):
    N_next[0] = 1
    N_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        N_next_Euler[i] = N[i] + dt*mu*(N[i]-(N[i])**2) #Euler
        N_next[i] = alpha[i] * N_next[i+1] + beta[i] # Crank-N
    N_next = N_next + N_next_Euler
#-implementar Euler explícito
#-Main, Condiciones iniciales, condiciones de bordes, busca y guarda las soluciones en cada paso, plotea

# Main

# setup
N_steps = 501
N_pasos_temporales = 400
gamma=0.001
#mu=1.5

h = 1 / (N_steps - 1)
# dt = h**2 / 2 # Este es el máximo teórico para el metodo explicito
dt = 0.01
r = (dt / 2 / h**2) * gamma

N = np.zeros(N_steps)
N_next = np.zeros(N_steps)
N_next_Euler = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_N(N, N_steps, h)


# Queremos guardar las soluciones en cada paso
N_solucion = np.zeros((N_pasos_temporales, N_steps))
N_solucion[0, :] = N.copy()

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(N, N_next, alpha, beta, N_steps)
    N = N_next.copy()
    N_solucion[i, :] = N.copy()


# Plots

# ejemplo 1

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 10):
    ax.plot(x, N_solucion[i, :])
ax.set_ylim(0, 1)

# ejemplo 2
# usar el plano x, t y plotear T en la 3a dimension
fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax2.pcolormesh(X, Y, N_solucion)

plt.show()
plt.draw()
