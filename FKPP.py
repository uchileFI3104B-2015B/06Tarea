# -*- coding: utf-8 -*-
"""
Created on Tue Nov 03 20:15:45 2015

@author: Administrador
"""

'''
Este script resuelve la ecuación de Fisher-KPP, usando
la discretización de Crank-Nicolson para la parte difusiva de la
ecuación y Euler explícito para la parte reactiva.
'''

import numpy as np
import matplotlib.pyplot as plt

def definir_condicion_inicial(n, h, L_x):
    for i in range(len(n)):
        n[i] = np.exp(-(h/L_x * i)**2 / 0.1)
    n[len(n)-1] = 0
    pass

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
