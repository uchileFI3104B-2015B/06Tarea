#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Este script contiene la la integracion de de La ecuación de Fisher-KPP
que describe el comportamiento animal mediante el metodo de euler y
cranck nicolson
'''


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

    def ecuacion_de_movimiento(n):
        '''
        Implementa la ecuación de movimiento
        '''
        fn = u * (n - n**2)
        return fn

    def sol_euler(y_actual,c):
        '''
        Toma la condición actual y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito.
        Devuelve self.y_actual y lo guarda en n_euler
        '''
        y_anterior = y_actual
        y_actual = y_anterior + dt * (ecuacion_de_movimiento(y_anterior))
        t_actual += dt
        n_euler[c] = y_actual
        return y_actual


def inicializa_n(n, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * h
        n[i] = np.exp( (- x**2)/0.1)
    n[0] = 1
    n[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * n[j+1] + (1-2*r) * n[j] + r * n[j-1]


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
N_steps = 41
N_pasos_temporales = 100

h = 1 / (N_steps - 1)
# dt = h**2 / 2 # Este es el máximo teórico para el metodo explicito
dt = 0.01
r = dt / 2 / h**2

y = 0.001
n = np.zeros(N_steps)
n_euler = np.zeros(N_pasos_temporales)
n_next = np.zeros(N_steps)
b = np.zeros(N_steps)
t0 = 1/y
x0 = 1

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
    n_next = n_next + sol_euler(n_euler[i], i)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()

for i in range (1, N_steps):
    Jupiter.avanza_euler(dt)
    resultados = Jupiter.y_actual
    x[i] = resultados[0]

# Plots

# ejemplo 1

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 10):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(0, 1)

# ejemplo 2
# usar el plano x, t y plotear T en la 3a dimension
fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax2.pcolormesh(X, Y, T_solucion)

plt.show()
plt.draw()
