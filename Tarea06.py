#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve un problema simple de diffusion en 1D.
La ecuación a resover es:
    dT/dt = d2T/dx2; T(0,x) = sin(pi * x); T(t, 0) = T(t, 1) = 0
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def inicializa_nu(T, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * h
        T[i] = np.exp(-10*x**2)
    T[0] = 1
    T[-1] = 0

def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * nu[ j + 1 ] + 2 * ( 1 - r ) * nu[ j ] + r * nu[ j-1 ]


def calcula_alpha_y_beta(alpha, beta, b, r, N_Steps):
    Aplus = - 1 * r
    Azero = ( 2 + 2 * r)
    Aminus = - 1 * r
    alpha[0] = 0
    beta[0] = 1
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Azero + Aminus * alpha[i - 1])
        beta[i] = (b[i] - Aminus*beta[i - 1]) / ( Aminus * alpha[i-1] + Azero)


def avanza_paso_crank(nu,nu_next, alpha, beta, N_steps):
    nu_next[0] = 1
    nu_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        nu_next[i] = alpha[i] * nu_next[i+1] + beta[i]

def avanza_paso_euler(nu, nu_next, dt, mu, N_steps):
    nu_next[0] = 1
    nu_next[-1] = 0
    for i in range(0,N_steps):
        nu_next[i] = nu[i] + mu * dt * ( nu[i] - nu[i]**2)



# setup
gamma=0.001
mu=1.5
N_steps = 500
N_pasos_temporales = 500
h = 1. / ( N_steps - 1 )
dt = 5./ ( N_pasos_temporales - 1 )
r = ( gamma * dt ) / ( 2 * h**2 )


nu = np.zeros(N_steps)
nu_euler = np.zeros(N_steps)
nu_crank = np.zeros(N_steps)
nu_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_nu(nu, N_steps, h)

# Queremos guardar las soluciones en cada paso
nu_solucion = np.zeros((N_pasos_temporales, N_steps))
nu_solucion[0, :] = nu.copy()

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    nu_crank=nu_solucion[i-1,:].copy()
    nu_euler=nu_solucion[i-1,:].copy()

    avanza_paso_crank(nu, nu_crank, alpha, beta, N_steps)

    avanza_paso_euler(nu, nu_euler, dt,mu, N_steps)

    nu = (nu_euler+nu_crank).copy()
    nu_solucion[i, :] = nu.copy()

# plots

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)


for i in range(0, N_pasos_temporales,200):
    ax.plot(x, nu_solucion[i, :])
#sax.set_ylim(0, 1)
plt.show()
plt.draw()


# ejemplo 2
# usar el plano x, t y plotear T en la 3a dimension
fig2 = plt.figure(2)
fig2.clf()
ax2 = fig2.add_subplot(111)
y = np.arange(0, N_pasos_temporales) * dt
X, Y = np.meshgrid(x, y)
ax2.pcolormesh(X, Y, nu_solucion)

plt.show()
plt.draw()
