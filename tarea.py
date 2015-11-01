#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pdb

'''Script para resolver numéricamente la ecuación de Fisher-KPP, usando
Crank Nicolson '''

e = 0.1     # paso temporal
h = 1 / 10 # paso espacial
GAMMA = 0.001
MU = 1.5
T_FIN = 4   # tiempo final
T_IN = 0    # tiempo inicial
X_FIN = 1   # extremo 2
X_IN = 0    # extremo 1
BORDE_1 = 1 # condicion de borde 1
BORDE_2 = 0 # condicion de borde 2
r = (GAMMA * e) / (2 * h ** 2)
A_menos = -r
A_0 = 1 + 2 * r
A_mas = -r

def funcion_inicial(x):
    '''función que representa la condicion inicial, recibe un x, y devuelve
    el valor de la función evaluada'''
    return np.exp(- x ** 2 / 0.1)


def iniciar_solucion(t, x):
    '''inicia la forma de la matriz con ceros, cada fila correspondera a una
    solución para n(x), donde la fila i es la solución para tiempo i.
    t es la cantidad de pasos temporales y x la de pasos espaciales'''
    return np.zeros((t+1, x+1))


def b(t, k):
    '''lado derecho (todo lo que depende del tiempo anterior
    al que se quiere calcular), para cada punto k, en el tiempo t'''
    b_valor = (solucion[t, k+1] * r + solucion[t, k-1] * r +
               solucion[t, k] * (MU * (1 - solucion[t, k]) - 2 * r))
    return b_valor


def poner_condiciones_iniciales():
    '''llena la fila 0 de la matriz con la condición inicial,
    retrona la nueva matriz'''
    x = np.linspace(X_IN, X_FIN, numero_pasos_x + 1)
    vfuncion_inicial = np.vectorize(funcion_inicial) # vectorizar funcion
    fila_inicial = vfuncion_inicial(x)
    solucion[0, :] = fila_inicial


def poner_condiciones_borde():
    '''define condiciones rígidas de borde para todo t'''
    solucion[:, 0] = BORDE_1
    solucion[:, numero_pasos_x] = BORDE_2


def encontrar_alfa_beta(t):
    ''' encuentra  todos los alfas y betas para el tiempo t dado'''
    # pdb.set_trace()
    N = numero_pasos_x - 1
    alfa = np.zeros(N + 1)         # iniciar parametros alfa y beta
    beta = np.zeros(N + 1)
    alfa[N] = 0                    # condiciones finales para alfa y beta
    beta[N] = BORDE_2

    for k in range(int(N), 0, -1): # iterar para llenar valores de alfa y beta
        alfa[k-1] = - A_menos / (A_mas * alfa[k] + A_0)
        beta[k-1] = (b(t, k) - A_mas * beta[k]) / (A_mas * alfa[k] + A_0)
    return alfa, beta


def iterar_crank_nicolson(t_actual):
    '''avanza un paso de la iteracion, es decir, llena la columna
    llamada con t_actual'''
    N = numero_pasos_x
    alfa, beta = encontrar_alfa_beta(t_actual)
    for k in range(0, int(N) - 1, int(N)):
        solucion[t_actual, k + 1] = alfa[k] * solucion [t_actual, k] + beta[k]


def mostar_resultado(sol):
    '''grafica el cada fila de sol vs x en el mismpo plot'''
    fig, ax = plt.subplots()
    x_values = np.linspace(X_IN, X_FIN, numero_pasos_x + 1)
    ax.plot(x_values, np.transpose(sol))
    plt.show()


# Main
numero_pasos_t = (T_FIN - T_IN) / e
numero_pasos_x = (X_FIN - X_IN) / h
solucion = iniciar_solucion(numero_pasos_t, numero_pasos_x)
poner_condiciones_iniciales()
poner_condiciones_borde()
iterar_crank_nicolson(1)
mostar_resultado(solucion)
