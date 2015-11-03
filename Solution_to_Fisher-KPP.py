#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Este script es para resolver la ecuación de Fisher-KPP
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

#-implementar Crank-N

def inicializa_N(N, N_steps, h):
    '''
    Rellena N con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean uno y cero.
    '''
    for i in range(N_steps):
        x = i * h
        N[i] = np.exp((-x**2) / 0.1)
    N[0] = 1
    N[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * N[j+1] + (1-2*r) * N[j] + r * N[j-1]


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde T(t, 0) = 1 ¿?
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(N, N_next, alpha, beta, N_steps):
    N_next[0] = 1
    N_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        N_next_Euler[i] = N[i] + dt*mu*(N[i]-(N[i])**2)
        N_next[i] = alpha[i] * N_next[i+1] + beta[i]

#-implementar Euler explícito
#-Main, Condiciones iniciales, condiciones de bordes, busca y guarda las soluciones en cada paso, plotea
