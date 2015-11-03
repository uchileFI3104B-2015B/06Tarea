#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Este script resuelve la ecuación de Fisher-Kolmogorov-Petrovsky-Piscounov,
usando la discretización de Crank-Nicolson para la parte difusiva de la
ecuación y Euler explícito para la parte reactiva.
'''

import numpy as np
import matplotlib.pyplot as plt


def set_condicion_inicial(n, h, Lx):
    for i in range(len(n)):
        n[i] = np.exp(-(h/Lx * i)**2 / 0.1)
    n[len(n)-1] = 0
    pass


def calculo_b(n, b, r, dt):
    '''
    Este método calcula el lado derecho de la solución usando Crank-Nicolson
    (es decir, los términos dependientes del instante temporal anterior).
    Recibe como parámetros los arreglos del campo n(x,t) en el insatante
    anterior y el arreglo b para guardar la información del lado derecho.
    El método modifica el arreglo y no retorna valor alguno.
    '''
    for i in range(1, len(b)-1):
        b[i] = r*n[i+1] + (1-2*r)*n[i] + r*n[i-1] + dt*(n[i] - n[i]**2)
    pass


def calc_alpha_beta(alpha, beta, b, r):
    '''
    Cálculo de arreglos para alpha y beta, que son el par de coeficientes que
    se usa para resolver cada paso temporal de la EDP.
    '''
    A_mas = -1 * r
    A_cer = (1 + 2 * r)
    A_men = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde n(t, 0) = 1
    for i in range(1, len(b)):
        alpha[i] = -A_mas / (A_cer + A_men*alpha[i-1])
        beta[i] = (b[i] - A_men*beta[i-1]) / (A_men*alpha[i-1] + A_cer)
    pass


def dt_step(n, n_sig, alpha, beta):
    '''
    Avanza un paso temporal en la solución.
    '''
    n_sig[0] = 1
    n_sig[-1] = 0
    for i in range(len(n) - 2, 0, -1):
        n_sig[i] = alpha[i] * n_sig[i+1] + beta[i]
    pass


mu = 1.5
gamma = 0.001
N_x = 501
# Cambio a sistema adimensional
L_x = 1.0/np.sqrt(gamma/mu)
h = 0.002*L_x
T = 7*mu
dt = (h**2)/4.0
N_t = int(np.ceil(T/dt))
r = dt/(2*h**2)

# Condición inicial y arreglos para guardar información
n = np.zeros(N_x)
set_condicion_inicial(n, h, L_x)
n_sig = np.zeros(N_x)
alpha = np.zeros(N_x)
beta = np.zeros(N_x)
b = np.zeros(N_x)
n_sol = np.zeros((N_t, N_x))
n_sol[0, :] = n.copy()

# Solución
for i in range(1, N_t):
    calculo_b(n, b, r, dt)
    calc_alpha_beta(alpha, beta, b, r)
    dt_step(n, n_sig, alpha, beta)
    n = n_sig.copy()
    n_sol[i, :] = n.copy()

# Plots
x = np.linspace(0, 1, N_x)
fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_t, 100):
    if i == 0:
        ax.plot(x, n_sol[i, :], color='b', label='$n(x,t_0)$ (tiempo fijo)')
    else:
        ax.plot(x, n_sol[i, :], color=(0, (1.0*i)/N_t, 1-(1.0*i)/N_t))
for j in range(50, N_x - 50, 50):
    ax.arrow(x[j], n_sol[100, j], 0.0,
             2.9*(n_sol[6900, j] - n_sol[100, j])/3 - 0.05, head_width=0.025,
             head_length=0.05, fc='r', ec='r')
ax.set_xlabel('$x$')
ax.set_ylabel('$n$')
ax.legend()
ax.set_title('Curvas de $n(x,t)$ para valores fijos de $t$')
ax.set_ylim([0, 1.15])
plt.savefig('FKPP.eps')
plt.show()
plt.draw()
