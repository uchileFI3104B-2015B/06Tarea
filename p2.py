#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Este script resuelve un problema simple de difusion en 1D.
La ecuación a resover es:
    dn/dt = gama * d2n/dx2 + un - un^3;
    C.I : n(0,x) = np.random.uniform(low=-0.3, high=0.3, size=steps_x);
    C.B : n(t, 0) = n(t, 1) = 0

La parte difusiva de la ecuación se utiliza el método de Crack-Nicolson
y la parte reactiva con el método de Euler explícito
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def inicializa_n(n, steps_x, h):
    '''
    Rellena n con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    func_random = np.random.uniform(low=-0.3, high=0.3, size=steps_x)
    for i in range(steps_x):
        n[i] = func_random[i]
    n[0] = 0
    n[-1] = 0


def calcula_b(b, steps_x, r, u=1.5):
    for j in range(1, steps_x - 1):
        # b[j] = r * n[j+1] + 2 * (1-r) * n[j] + r * n[j-1]
        b[j] = (r * n[j+1] + (1-2*r) * n[j] + r * n[j-1] +
                n[j] * (dt * u - dt * u * n[j]**2))


def calcula_alfa_y_beta(alfa, beta, b, r, steps_x):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alfa[0] = 0
    beta[0] = 0  # viene de la condicion de borde n(t, 0) = 1
    for i in range(1, steps_x):
        alfa[i] = -Aplus / (Acero + Aminus*alfa[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alfa[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alfa, beta, steps_x, steps_t):
    n_next[0] = 0
    n_next[-1] = 0

    for i in range(steps_x - 2, 0, -1):
        n_next[i] = alfa[i] * n_next[i+1] + beta[i]


# Main

# setup
seteo = input('Valor semilla?: ')
np.random.seed(seteo)
gamma = 0.001
mu = 1.5
ti = 0
tf = 20
xi = 0
xf = 1

steps_t = 1000
steps_x = 500

dt = (tf - ti) / (steps_t - 1)
h = (xf - xi) / (steps_x - 1)

# dt = h**2 / 2 # Este es el máximo teórico para el metodo explicito
# Máximo teórico no funciona

r = dt / (2 * h**2)
r = r * gamma

n = np.zeros(steps_x)
n_next = np.zeros(steps_x)

b = np.zeros(steps_x)
alfa = np.zeros(steps_x)
beta = np.zeros(steps_x)

inicializa_n(n, steps_x, h)

# Queremos guardar las soluciones en cada paso
n_solucion = np.zeros((steps_t, steps_x))
n_solucion[0, :] = n.copy()

for i in range(1, steps_t):
    calcula_b(b, steps_x, r)
    calcula_alfa_y_beta(alfa, beta, b, r, steps_x)
    avanza_paso_temporal(n, n_next, alfa, beta, steps_x, steps_t)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()


# Plot
x = np.linspace(0, xf, steps_x)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, steps_t, tf):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(-1, 1)

plt.xlabel('x', fontsize=14)
plt.ylabel('n (x)', fontsize=14)
plt.title("Tiempo final = {tf}, Semilla = {seteo}".format(tf=tf, seteo=seteo),
          fontsize=16)
plt.savefig("p2_tf{tf}_setrandom{seteo}.png".format(tf=tf, seteo=seteo))

plt.show()
plt.draw()
