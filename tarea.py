#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pdb

'''Script para resolver numéricamente la ecuación de Fisher-KPP, usando
Crank Nicolson '''

def funcion_inicial(x):
    '''función que representa la condicion inicial, recibe un x, y devuelve
    el valor de la función evaluada'''
    return np.exp(- x ** 2 / 0.1)


def poner_condiciones_iniciales(S, N, h):
    '''llena la fila 0 de la matriz con la condición inicial,
    retrona la nueva matriz'''
    x = np.linspace(X_IN, X_FIN, N)
    vfuncion_inicial = np.vectorize(funcion_inicial) # vectorizar funcion
    fila_inicial = vfuncion_inicial(x)
    S = fila_inicial
    S[0] = BORDE_1
    S[-1] = BORDE_2
    return S


def calcula_b(b, N, r):
    '''lado derecho (todo lo que depende del tiempo anterior
    al que se quiere calcular), para cada punto k'''
    for k in range(1, N - 1):
        b[k] = (S[k+1] * r + S[k-1] * r +
               S[k] * (e * MU * (1 - S[k]) + 1 - 2 * r))


def encontrar_alfa_beta(alfa, beta, b, r,N ):
    ''' encuentra  todos los alfas y betas para el tiempo t dado'''
    A_mas = -1 * r
    A_menos = -1* r
    A_0 = (1 + 2 * r)
    alfa[0] = 0                # condiciones finales para alfa y beta
    beta[0] = BORDE_2

    for k in range(1, N): # iterar para llenar valores de alfa y beta
        alfa[k] = -1 * A_mas / (A_menos * alfa[k-1] + A_0)
        beta[k] = (b[k] - A_menos * beta[k-1]) / (A_menos * alfa[k-1] + A_0)


def avance_temporal(S, S_next, alfa, beta, N):
    '''avanza un paso de la iteracion temporal, es decir, llena la columna
    llamada con t'''
    S_next[0] = BORDE_1
    S_next[-1] = BORDE_2
    for k in range(int(N) - 2, 0, -1):
        S_next[k] = alfa[k] * S_next[k+1] + beta[k]


def mostar_resultado(sol):
    '''grafica el cada fila de sol vs x en el mismpo plot'''
    fig, ax = plt.subplots()
    x_values = np.linspace(X_IN, X_FIN, numero_x)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)
    for i in range(0, numero_t):
        ax.plot(x_values, sol[i, :], color="r")
    plt.show()

# Main
# Setup
T_FIN = 4       # tiempo final
T_IN = 0        # tiempo inicial
X_FIN = 1       # extremo 2
X_IN = 0        # extremo 1
numero_t = 100       # cantidad de valores en t
numero_x = 500   # cantidad de valores en x
numero_pasos_t = numero_t - 1
numero_pasos_x = numero_x - 1
e = (T_FIN - T_IN) / (numero_pasos_t)
h = (X_FIN - X_IN) / (numero_pasos_x)
GAMMA = 0.001
MU = 1.5
BORDE_1 = 1     # condicion de borde 1
BORDE_2 = 0     # condicion de borde 2
r = (GAMMA * e) / (2 * h ** 2)

S = np.zeros(numero_x)
S_next = np.zeros(numero_x)

b = np.zeros(numero_x)
alfa = np.zeros(numero_x)
beta = np.zeros(numero_x)

S = poner_condiciones_iniciales(S, numero_x, h)

S_solucion = np.zeros((numero_t, numero_x))
S_solucion[0, :] = S.copy()

for i in range(1, numero_t):
    calcula_b(b, numero_x, r)
    encontrar_alfa_beta(alfa, beta, b, r, numero_x)
    avance_temporal(S, S_next, alfa, beta, numero_x)
    S = S_next.copy()
    S_solucion[i, :] = S.copy()

mostar_resultado(S_solucion)
