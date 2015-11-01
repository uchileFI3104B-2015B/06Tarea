#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def set_condicion_inicial(n, h, Lx):
    for i in range(len(n)):
        n[i] = np.exp(-(h/Lx * i)**2 / 0.1)
    pass


mu = 1.5
gamma = 0.001
N_x = 501
# Cambio a sistema adimensional
L_x = 1.0/np.sqrt(gamma/mu)
h = 0.002*L_x
T = 10*mu
dt = (h**2)/4.0

n_in = np.zeros(N_x)
set_condicion_inicial(n_in, h, L_x)
