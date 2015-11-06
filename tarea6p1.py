'''
Este script resuelve la ecuación de Fisher-KPP mediante
el método de Crank-Nicolson.
La ecuación a resolver es:
dn/dt = gamma * d2n/dx2 + mu* n - mu* n^2
Condiciones de borde e incial:
n(t, 0) = 1 ; n(t, 1) = 0 ; n(0, x) = exp(-x^2 / 0.1)
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def inicializa_n(n, N_steps, dx): #TO DO: definir N_steps, dx y n
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * dx
        n[i] = np.exp((-x**2) / 0.1)
    n[0] = 1
    n[-1] = 0
    pass

def calcular_b(b, N_steps, r):
    '''
    Calcula la parte derecha de la ecuación, que depende
    del instante anterior.
    '''
    for j in range(1, N_steps - 1):
        b[j] = r*n[j+1] + (1 - 2*r)*n[j] + r*n[j-1] + dt*mu*(n[j] - n[j]**2) #TO DO: definir r, mu y dt
    pass

def alpha_beta(alpha, beta, b, r, N_steps): #TO DO: definir alpha y beta
    '''
    Calcula los valores de alpha y beta.
    '''
    A_plus = -1 * r
    A_zero = (1 + 2 * r)
    A_minus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde n(t, 0) = 1
    for i in range(1, N_steps):
        alpha[i] = -A_plus / (A_zero + A_minus * alpha[i-1])
        beta[i] = (b[i] - A_minus * beta[i-1]) / (A_minus * alpha[i-1] + A_cero)
    pass

def time_step(n, n_next, alpha, beta): #TO DO: definir alpha, beta y n_next
    '''
    Avanza un instante en el tiempo.
    '''
    n_next[0] = 1 #condición de borde n(t,0) = 1
    n_next[-1] = 0 #condición de borde n(t,0) = 0
    for i in range(N_steps - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]
    pass
