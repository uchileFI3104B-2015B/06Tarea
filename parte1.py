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
