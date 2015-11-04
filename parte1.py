#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve la ecuacion de Fisher-KPP
La ecuación a resover es:
    dn/dt = gamma * d2n/dx2 + mu * n - mu * n ** 2
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# funciones base


def cond_ini(x):
    '''
    establece las condiciones iniciales del problema
    '''
    pass


def calc_b(b, M, r):
    '''
    calcula el valor de b para el metodo de crank nicolson
    '''
    for j in range(1, M - 1):
        b[j] = r * n[]
    pass


def calc_alpha_beta(alpha, beta, b, r, M):
    '''
    obtiene los valores de alpha y beta para el metodo de crank nicolson
    '''
    pass


def avanza_tiempo(n, n_sig, alpha, beta, M):
    '''
    obtiene el comportamiento de la funcion al tiempo siguiente
    '''
    pass



# main
# inicializacion

# iteracion

# plots
