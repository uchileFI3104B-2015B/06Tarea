#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def set_condicion_inicial(n, h, Lx):
    for i in range(len(n)):
        n[i] = np.exp(-(h/Lx * i)**2 / 0.1)
    pass

def calculo_b(n, b, r, dt):
    for i in range(1, len(b)-1):
        b[i] = r*n[i+1] + (1-2*r)*n[i] + r*n[i-1] + dt*(n[i] + n[i]**2)

def calc_alpha_beta():
    pass

def dt_step():
    pass

mu = 1.5
gamma = 0.001
N_x = 501
# Cambio a sistema adimensional
L_x = 1.0/np.sqrt(gamma/mu)
h = 0.002*L_x
T = 10*mu
dt = (h**2)/4.0
N_t = np.ceil(T/dt)

r = dt/(2*h**2)
n = np.zeros(N_x)
set_condicion_inicial(n, h, L_x)
n_sig = np.zeros(N_x)
alpha = np.zeros(N_x)
beta = np.zeros(N_x)
b = np.zeros(N_x)
n_sol = np.zeros((N_t,N_x))
n_sol[0, :] = n.copy()

for i in range(1, N_t):
    calculo_b(n, b, r, dt)
    calc_alpha_beta(alpha, beta, b, r, N_x)
    dt_step(n, n_sig, alpha, beta, N_x)
    n = n_sig.copy()
    n_sol[i, :] = n.copy()
