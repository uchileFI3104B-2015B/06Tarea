#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Este script es para resolver la ecuación de Newell-Whitehead-Segel
'''

from __future__ import division
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

np.random.seed(10)
#-implementar Crank-N

def inicializa_N(N, N_steps, h):
    '''
    Rellena N con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean uno y cero.
    '''
    N = np.random.uniform(low=-0.3, high=0.3, size=N_steps)
    print N[0], N[1],N[2]
    N[0] = 1
    N[-1] = 0
    return N

def calcula_b(b, N_steps, r, mu=1.5):
    for j in range(1, N_steps - 1):
        b[j] = r * N[j+1] + (1-2*r) * N[j] + r * N[j-1] + dt*mu*(N[j]-(N[j])**3)


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
        N_next[i] =alpha[i] * N_next[i+1] + beta[i]  # Crank-N

# Main
# setup
N_steps = 501
N_pasos_temporales = 1000
gamma=0.001
#mu=1.5

h = 1 / (N_steps - 1)
dt = 0.01
r = (dt / 2 / h**2) * gamma

N = np.zeros(N_steps)
N_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

N = inicializa_N(N, N_steps, h)
print N[0],N[1],N[2]
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

# PLOT 1

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 10):
    ax.plot(x, N_solucion[i, :],'r')
#ax.set_ylim(0, 1)

plt.xlabel('Posicion')
plt.ylabel('Densidad de la especie')
plt.title("Densidad de la especie en funcion de la posicion entre t = 0 y t = 4")
plt.savefig("Imagen_1")
#Plot 2

fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
from matplotlib.collections import LineCollection
line_segments = LineCollection([list(zip(x, ys)) for ys in N_solucion],
                               linewidths=(0.5, 1, 1.5, 2),
                               linestyles='solid')
line_segments.set_array(x)
ax2.add_collection(line_segments)
fig2 = plt.gcf()
plt.sci(line_segments)

plt.xlabel('Posicion')
plt.ylabel('Densidad de la especie')
plt.title("Densidad de la especie en funcion de la posicion entre t = 0 y t = 4")
plt.savefig("Imagen_2")
# PLOT 3
# usar el plano x, t y plotear N en la 3a dimension

fig3 = plt.figure(3)
fig3.clf()
ax3 = fig3.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax3.pcolormesh(X, Y, N_solucion)


plt.xlabel('Posicion')
plt.ylabel('Tiempo')
plt.title("Plano distancia-tiempo con densidad en la tercera dimension")
plt.savefig("Imagen_3")
plt.show()
plt.draw()
