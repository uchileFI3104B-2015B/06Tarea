#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt



'''Script para resolver numéricamente la ecuación de Newell-Whitehead-Segel, usando
Crank Nicolson.
dn/dt = gamma * d2n/dx2 + mu* n - mu * n^3 '''

np.random.seed(10)

def poner_condiciones_iniciales(S, N, h):
    '''llena la fila 0 de la matriz con la condición inicial,
    retrona la nueva matriz'''
    S = np.random.uniform(low=-0.3, high=0.3, size=N)
    S[0] = CdB1
    S[-1] = CdB2
    return S


def calcula_b(b, N, r):
    '''lado derecho (todo lo que depende del tiempo anterior
    al que se quiere calcular), para cada punto k'''
    for i in range(1, N - 1):
        b[i] = (S[i+1] * r + S[i-1] * r +
                S[i] * (eps * Mu * (1 - S[i]**2) + 1 - 2 * r))


def encontrar_alfa_beta(alfa, beta, b, r,N ):
    ''' encuentra  todos los alfas y betas para el tiempo t dado'''
    A_mas = -1 * r
    A_menos = -1* r
    A_0 = (1 + 2 * r)
    alfa[0] = 0
    beta[0] = CdB1

    for i in range(1, N):
        alfa[i] = -1 * A_mas / (A_menos * alfa[i-1] + A_0)
        beta[i] = (b[i] - A_menos * beta[i-1]) / (A_menos * alfa[i-1] + A_0)


def avance_temporal(S, S_next, alfa, beta, N):
    '''avanza un paso de la iteracion temporal, es decir, llena la columna
    llamada con t'''
    S_next[0] = CdB1
    S_next[-1] = CdB2
    for i in range(int(N) - 2, 0, -1):
        S_next[i] = alfa[i] * S_next[i+1] + beta[i]





Tf = 10
To = 0
Xf= 1
Xo = 0
num_t = 100
num_x = 500
num_pasos_t = num_t - 1
num_pasos_x = num_x - 1
eps = (Tf - To) / (num_pasos_t)
h = (Xf - Xo) / (num_pasos_x)
Gamma = 0.001
Mu = 1.5
CdB1 = 0
CdB2 = 0
r = (Gamma * eps) / (2 * h ** 2)

S = np.zeros(num_x)
S_next = np.zeros(num_x)

b = np.zeros(num_x)
alfa = np.zeros(num_x)
beta = np.zeros(num_x)

S = poner_condiciones_iniciales(S, num_x, h)

Ssol = np.zeros((num_t, num_x))
Ssol[0, :] = S.copy()

for i in range(1, num_t):
    calcula_b(b, num_x, r)
    encontrar_alfa_beta(alfa, beta, b, r, num_x)
    avance_temporal(S, S_next, alfa, beta, num_x)
    S = S_next.copy()
    Ssol[i, :] = S.copy()




fig, ax = plt.subplots()
val_x = np.linspace(Xo, Xf, num_x)
ax.set_xlabel("X en unidades espaciales")
ax.set_ylabel("N (X,T)")
ax.set_title("N (X,T) para tiempos entre t=0 y t=10")
for i in range(0, num_t):
    ax.plot(val_x, Ssol[i, :], color="green")


plt.show()
plt.draw()
