# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

'''
Este script resuelve la ecuación de Fisher-KPP, usando
la discretización de Crank-Nicolson para la parte difusiva de la
ecuación y Euler explícito para la parte reactiva.
'''


def definir_condiciones_borde(n, N, h):
    '''llena la fila 0 de la matriz con la condición inicial,
    retrona la nueva matriz'''
    for i in range(N):
        x = i * h
        n[i] = np.exp(-1 * (x ** 2) / 0.1)
    n[0] = B1
    n[-1] = B2
    return n

def calcular_b(b, n, r, N, mu):
    '''Calcula el lado derecho (todo lo que depende del tiempo anterior
    al que se quiere calcular), para cada punto i'''
    for i in range(1, N-1):
          b[i] = (r * n[i+1] + r * n[i-1] + 
                  n[i] * (dt * mu * (1 - n[i]) + 1 - 2 * r))

def calcular_alpha_y_beta(alpha, beta, b, r, N):
    '''Calcula todos los alphas y betas para el tiempo t dado'''
    Amas = -1 * r
    Acero = (1 + 2 * r)
    Amenos = -1 * r
    alpha[0] = 0
    beta[0] = B1  # debido a que condicion de borde n(t, 0) = 1
    for i in range(1, N):
        alpha[i] = -1 * Amas / (Acero + Amenos * alpha[i-1])
        beta[i] = (b[i] - Amenos * beta[i-1]) / (Amenos * alpha[i-1] + Acero)
        
def avance_tiempo(n_t, n_tsig, alpha, beta, N):
    '''
    Avanza un paso temporal en la solución.
    '''
    n_tsig[0] = B1 #Condiciones de borde
    n_tsig[-1] = B2
    for i in range(int(N) - 2, 0, -1):
        n_tsig[i] = alpha[i] * n_tsig[i+1] + beta[i]

#Main

#Inicializacion
B1=1
B2=0
mu = 1.5
gamma = 0.001
N_x = 501
L_x = 1
h = 1.0/500.0
T = 10
dt = 0.1
N_t = 101
r = (gamma * dt) / (2 * h ** 2)
#C.Iniciales y definicion de arreglos
n_t = np.zeros(N_x)

n_t = definir_condiciones_borde(n_t, N_x, h)
n_tsig = np.zeros(N_x)

alpha = np.zeros(N_x)
beta = np.zeros(N_x)
b = np.zeros(N_x)

n_sol = np.zeros((N_t, N_x))
n_sol[0, :] = n_t.copy()
#Iteración
for k in range(1, N_t):
    calcular_b(b, n_t, r, N_x, mu)
    calcular_alpha_y_beta(alpha, beta, b, r, N_x)
    avance_tiempo(n_t, n_tsig, alpha, beta, N_x)
    n_t = n_tsig.copy()
    n_sol[k, :] = n_t.copy()
    
#Graficar
fig = plt.figure(1, figsize=(7,6))
fig.clf()
ax = fig.add_subplot(111)
x_values = np.linspace(0, 1, N_x)
ax.set_ylim(0, 1)
ax.set_xlim(0, 1)
ax.set_xlabel("X en unidades arbitrarias de espacio")
ax.set_ylabel("Densidad en unidades arbitrarias")
ax.set_title("Densidad en el espacio para tiempo entre t=0 y t=10")
for k in range(0, N_t, 10):
    ax.plot(x_values, n_sol[k, :], label="t="+str(k*dt))
plt.legend(loc='lower left')
fig.savefig("FKPP.jpg")
plt.show()