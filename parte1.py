'''
Este script resuelve un problema simple de diffusion en 1D.
La ecuación a resover es:
    dn/dt = gama * d2n/dx2 + un - un^2; n(0,x) = np.exp(-x**2/0.1); n(t, 0) = 1 n(t, 1) = 0
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def inicializa_T(n, N_steps, h):
    '''
    Rellena n con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * h
        n[i] = np.exp(-x**2/0.1)
    n[0] = 1  # c de borde
    n[-1] = 0  # condiciones de borde


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = (r * n[j+1] + (2-2*r) * n[j] + r * n[j-1] +
                dt * mu * n[j] - dt * (mu) * n[j] * n[j])  # solucion


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (2 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):
    T_next[0] = 1
    T_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]
    # T = T_next.copy() # Esto no funciona, hacerlo fuera de la funcion


# Main

# setup
gamma = 0.001
mu = 1.5
t_inicial = 0
t_final = 45
x_inicial = 0
x_final = 1

N_steps = 500
N_pasos_temporales = 1000

dt = (t_final - t_inicial) / (N_pasos_temporales - 1)
h = (x_final - x_inicial) / (N_steps-1)
r = (dt / (2 * h**2))*gamma
