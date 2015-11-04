# -*- coding: utf-8 -*-
"""
Created on Tue Nov 03 20:15:45 2015

@author: Administrador
"""

'''
Este script resuelve la ecuación de Fisher-iPP, usando
la discretización de Crani-Nicolson para la parte difusiva de la
ecuación y Euler explícito para la parte reactiva.
'''

import numpy as np
import matplotlib.pyplot as plt

def definir_condiciones_borde(n, h, L_x):
    for i in range(len(n)):
        n[i] = np.exp(-(h/L_x * i)**2 / 0.1)
    n[len(n)-1] = 0
    pass

def calcular_b(n, b, r, dt):
    '''Calcula el lado derecho (todo lo que depende del tiempo anterior
    al que se quiere calcular), para cada punto i'''
    for i in range(1, len(b)-1):
        b[i] = r*n[i+1] + (1-2*r)*n[i] + r*n[i-1] + dt*(n[i] - n[i]**2)
    pass

def calcular_alpha_y_beta(alpha, beta, b, r):
    '''Calcula todos los alphas y betas para el tiempo t dado'''
    Amas = -1 * r
    Acero = (1 + 2 * r)
    Amenos = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde n(t, 0) = 1
    for i in range(1, len(b)):
        alpha[i] = -Amas / (Acero + Amenos*alpha[i-1])
        beta[i] = (b[i] - Amenos*beta[i-1]) / (Amenos*alpha[i-1] + Acero)

#Main

#Inicializacion
mu = 1.5
gamma = 0.001
N_x = 501
#Adimensionalizacion
L_x = 1.0/np.sqrt(gamma/mu)
h = 0.002 * L_x
T = 7 * mu
dt = (h**2) / 4.0
N_t = int(np.ceil(T/dt))
r = dt / (2*h**2)
#C.Iniciales y definicion de arreglos
n_t = np.zeros(N_x)
definir_condiciones_borde(n_t, h, L_x)
n_tsig = np.zeros(N_x)
alpha = np.zeros(N_x)
beta = np.zeros(N_x)
b = np.zeros(N_x)
n_sol = np.zeros((N_t, N_x))
n_sol[0, :] = n_t.copy()
