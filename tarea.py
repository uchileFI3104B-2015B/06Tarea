#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import pdb

'''Script para resolver numéricamente la ecuación de Fisher-KPP, usando
Crank Nicolson '''
e = 0.1     # paso temporal
h = 1 / 10 # paso espacial
GAMMA = 0.001
NU = 1.5
T_FIN = 4   # tiempo final
T_IN = 0    # tiempo inicial
X_FIN = 1   # extremo 2
X_IN = 0    # extremo 1
BORDE_1 = 1 # condicion de borde 1
BORDE_2 = 0 # condicion de borde 2

def funcion_inicial(x):
    '''función que representa la condicion inicial, recibe un x, y devuelve
    el valor de la función evaluada'''
    return np.exp(- x ** 2 / 0.1)



def iniciar_solucion(t, x):
    '''inicia la forma de la matriz con ceros, cada fila correspondera a una solución
    para n(x), donde la fila i es la solución para tiempo i'''
    return np.zeros((t+1, x+1))


def poner_condiciones_iniciales():
    '''llena la fila 0 de la matriz con la condición inicial,
    retrona la nueva matriz'''
    x = np.linspace(X_IN, X_FIN, numero_pasos_x + 1)
    vfuncion_inicial = np.vectorize(funcion_inicial) # vectorizar funcion
    fila_inicial = vfuncion_inicial(x)
    solucion[0, :] = fila_inicial


def iterar_crank_nicolson():
    '''avanza un paso de la iteracion, es decir, llena la columna siguiente'''


def mostar_resultado():
    '''  grafica el resultado'''


# Main
numero_pasos_t = (T_FIN - T_IN) / e
numero_pasos_x = (X_FIN - X_IN) / h
solucion = iniciar_solucion(numero_pasos_t, numero_pasos_x)
poner_condiciones_iniciales()
# solucion = iterar_crank_nicolson()
mostar_resultado()
