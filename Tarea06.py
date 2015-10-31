#!/usr/bin/env python
# -*- coding: utf-8 -*-


'''
Este script resuelve un problema simple de diffusion en 1D.
La ecuación a resover es:
    dT/dt = d2T/dx2; T(0,x) = sin(pi * x); T(t, 0) = T(t, 1) = 0
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * T[j+1] + 2*(1-r) * T[j] + r * T[j-1]


def calcula_alpha_y_beta(alpha, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (2 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_crank(nu_next, alpha, beta, N_steps):
    nu_next[0] = 1
    nu_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        nu_next[i] = alpha[i] * nu_next[i+1] + beta[i]

def avanza_paso_euler(nu,nu_next,h,mu, N_steps):
    nu_next[0] = 1
    nu_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        nu_next[i] = nu[i] + mu * h * ( nu_next[i] - nu_next[i]**2)
