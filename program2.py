#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Se resuelve la ecuacion de Fisher-KPP. La parte difusiva de la ecuación
se calcula con el metodo de Crank Nicolson y la parte reactiva
con el método de Euler Explicito
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(58)


def inicializa_n(ci, n_steps, h):
    ran = np.random.uniform(low=-0.3, high=0.3, size=n_steps)
    for i in range(n_steps):
        n[i] = ran[i]
    ci[0] = 0
    ci[-1] = 0


def calcula_b(b, n_steps, r, u=1.5):
    for j in range(1, n_steps - 1):
        b[j] = r * n[j+1] + (1 - 2 * r) * n[j] + r * n[j-1] + n[j] * (dt * u -
                        dt * u * n[j] ** 2)


def calcula_alfa_y_beta(alfa, beta, b, r, n_steps):
    Aplus = -r
    Acero = (1 + 2 * r)  # acero xd
    Aminus = -r
    alfa[0] = 0
    beta[0] = 0
    for i in range(1, n_steps):
        alfa[i] = - Aplus / (Acero + Aminus * alfa[i-1])
        beta[i] = (b[i] - Aminus * beta[i-1])/(Aminus * alfa[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alfa, beta, n_steps):
    n_next[0] = 0    # pq n(t, 0) = 0
    n_next[-1] = 0   # pq n(t, 1) = 0
    for i in range(n_steps-2, 0, -1):
        n_next[i] = alfa[i] * n_next[i+1] + beta[i]

# Ploteo


def graficar(x, t, sol):
    x = np.linspace(0, 1, x)

    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)

    for i in range(0, t, 10):
        ax.plot(x, n_solucion[i, :])
    ax.set_ylim(-1, 1)

    plt.xlabel('x')
    plt.ylabel('n')
    plt.title('Grafica de n(x,t) para randomseed=58')

    fig.savefig('n2(x,t)_2.png')

    plt.show()
    plt.draw()

# condiciones del problema

ti, tf = 0, 10  # ahora son 10 sec
xi, xf = 0, 1

steps_t = 1000
steps_x = 500

dt = (tf - ti) / (steps_t - 1)
dx = (xf - xi) / (steps_x - 1)

r = dt / (2 * dx ** 2)
r = r * 0.001  # gamma es 0.001

n = np.zeros(steps_x)
n_next = np.zeros(steps_x)
b = np.zeros(steps_x)
alfa = np.zeros(steps_x)
beta = np.zeros(steps_x)

inicializa_n(n, steps_x, dx)

# Queremos guardar las soluciones en cada paso
n_solucion = np.zeros((steps_t, steps_x))
n_solucion[0, :] = n.copy()


for i in range(1, steps_t):
    calcula_b(b, steps_x, r)
    calcula_alfa_y_beta(alfa, beta, b, r, steps_x)
    avanza_paso_temporal(n, n_next, alfa, beta, steps_x)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()


graficar(steps_x, steps_t, n_solucion)
