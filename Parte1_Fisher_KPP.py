#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve la ecuacion de difusion-reaccion Fisher-KPP
dn/dt=gamma d2n/dx2 + mu*n-mu*n2
Utiliza metodo de Crank-Nicolson para la parte de difusion de la EDP
y euler explicito para la parte de reaccion
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def inicializa_n(n, N_pasos, h):
    '''
    Rellena n con las condiciones iniciales y de bordes del problema
    '''
    for i in range(N_pasos):
        x = i * h
        n[i] = np.exp(- x**2/0.1)
    n[0] = 1
    n[-1] = 0


def calcula_b(b, N_pasos, r):
    """
    calcula termino b de acuerdo a la discretizacion
    """
    for j in range(1, N_pasos - 1):
        b[j] = r * n[j+1] + (1-2*r) * n[j] + r * n[j-1] + dt*(n[j] - n[j]**2)


def calcula_alpha_y_beta(alhpa, beta, b, r, N_pasos):
    """
    calcula terminos alpha y beta de la discretizacion
    """
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde n(t, 0) = 1
    for i in range(1, N_pasos):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alpha, beta, N_pasos):
    """
    calcula solución para el siguiente paso de tiempo
    """
    n_next[0] = 1
    n_next[-1] = 0
    for i in range(N_pasos - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]

# Estructura

# parametros
gamma = 0.001
mu = 1.5

# setup
N_pasos = 501
N_pasos_temporales = 100
h = 1 / (N_pasos - 1)
T_total = 15.0
dt = T_total/N_pasos_temporales
r = (gamma)*dt / 2 / h**2
n = np.zeros(N_pasos)
n_next = np.zeros(N_pasos)
b = np.zeros(N_pasos)
alpha = np.zeros(N_pasos)
beta = np.zeros(N_pasos)

inicializa_n(n, N_pasos, h)


# Se guardan las soluciones en cada paso
n_solucion = np.zeros((N_pasos_temporales, N_pasos))
n_solucion[0, :] = n.copy()

# se integra la ecuacion

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_pasos, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_pasos)
    avanza_paso_temporal(n, n_next, alpha, beta, N_pasos)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()

# Plots

# Grafico 2D densidad vs posicion para en distintos tiempos

x = np.linspace(0, 1, N_pasos)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 1):
    ax.plot(x, n_solucion[i, :], color='b')
plt.xlabel('posicion x')
plt.ylabel('densidad n(x,t)')
plt.title('densidad vs posicion para 0 < t < 15 ')

# Grafico 3D densidad vs posicion vs tiempo

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax2.plot_wireframe(X, Y, n_solucion, rstride=5, cstride=5)
ax2.set_title('densidad en el espacio para 0 < t < 15')
ax2.set_xlabel('posicion x')
ax2.set_ylabel('tiempo t')
ax2.set_zlabel('densidad n(x,t)')
plt.show()
plt.draw()
