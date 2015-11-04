##############################################################################
'''
Metodos Numericos para la Ciencia e Ingenieria
FI3104-1
Tarea 6
Maximiliano Dirk Vega Aguilera
18.451.231-9
'''
##############################################################################
##############################################################################

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

##############################################################################
##############################################################################
#funciones

def inicializa_n(n, N_pasos, h):
    '''
    Rellena n con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_pasos):
        x = i * h
        n[i] = np.exp((-x**2)/0.1)
    n[0] = 1
    n[1] = 0


def calcula_b(b, N_pasos, r):
    for j in range(1, N_pasos - 1):
        b[j] = r * n[j+1] + (1-2*r) * n[j] + r * n[j-1]


def calcula_alpha_y_beta(alhpa, beta, b, r, N_pasos):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_pasos):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alpha, beta, N_pasos):
    n_next[0] = 0
    n_next[-1] = 0
    for i in range(N_pasos - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]

##############################################################################
##############################################################################
#Desarrollo P1
'''
dn/dt = g*d2n/dx2 + mn - mn2   n(x,t)
crank-nicolson para el gradiente, euler explicito para mu
0<x<1
gamma = 0.001
mu = 1.5
discretizar el espacio en unos 500 pasos
condiciones:
n(t,0) = 1
n(t,1) = 0
n(0,x) = exp(-x**2 / 0.1)
integrar para almenos t=4
'''
N_pasos = 500   #numero de pasos espaciales
h = 1./N_pasos  #tamanho del paso espacial
e = 0.005       #tamanho del paso espacial




##############################################################################
##############################################################################
#Desarrollo P2
'''
dn/dt = g*d2n/dx2 + mn - mn3   n(x,t)
'''

##############################################################################
##############################################################################
#Graficos

##############################################################################
##############################################################################
#zona de pruebas







##############################################################################
