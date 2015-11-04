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


def cond_ini(n, M, h):
    '''
    establece las condiciones iniciales del problema
    '''
    for i in range(N):
        x = i * h
        n[i] = np.exp(- (x ** 2) / 0.1)
    pass


def calc_b(b, M, r, h):
    '''
    calcula el valor de b para el metodo de crank nicolson
    '''
    for j in range(1, M - 1):
        b[j] = r * n[j+1] + (1 - 2 * r) * n[j] +  r * n[j-1] +
               h * n[i] * (1 - n[i])
    pass


def calc_alpha_beta(alpha, beta, b, r, M):
    '''
    obtiene los valores de alpha y beta para el metodo de crank nicolson
    '''
    Amas = - 1 * r
    Acero = (1 + 2 * r)
    Amenos = - 1 * r
    alpha[0] = 0
    beta[0] = 1  # condiciones iniciales
    for i in range(1, M):
        alpha[i] = - Amenos / (Acero + Amenos * alpha[i-1])
        beta[i] = (b[i] - Amenos * beta[i-1]) / (Amenos * alpha[i-1] + Acero)
    pass


def avanza_tiempo(n, n_sig, alpha, beta, M):
    '''
    obtiene el comportamiento de la funcion al tiempo siguiente
    '''
    n_sig[0] = 1  # condiciones iniciales
    n_sig[-1] = 0
    for i in range(M - 2, 0, -1):
        n_sig[i] = alpha[i] * n_sig[i+1] + beta[i]
    pass



# main
# inicializacion

# iteracion

# plots
