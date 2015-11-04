
# coding: utf-8

# In[31]:

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# P1 Tarea 6


class Crank_Nicolson(object):
    def __init__(self, condiciones_borde, condiciones_iniciales):
        '''Inicializa la clase con los valores
        de las condiciones de borde, condiciones
        iniciales, mu, gamma, etc.'''
        self.cb = condiciones_borde  # n(x=0,t)=1 y n(x=1,t)=0
        self.ci = condiciones_iniciales  # n(x,0)=exp(-x**2/0.1)
        self.mu = 1.5
        self.gamma = 0.001
        self.t_actual = 0
        self.numinterx = 499
        self.numptosx = 500
        self.dx = 1/self.numinterx

    # matrices tridiagonales

    def calcula_d(self, n, dt, r, pregunta):
        d = np.zeros(self.numptosx)
        if pregunta == 1:
            for i in range(1, self.numptosx - 1):
                d[i] = r * n[i+1] + (1-2*r) * n[i] +                        r * n[i-1] + dt*(n[i] - n[i]**2)
        else:
            for i in range(1, self.numptosx - 1):
                d[i] = r * n[i+1] + (1-2*r) * n[i] +                        r * n[i-1] + dt*(n[i] - n[i]**3)
        return d

    def calcula_alpha_y_beta(self, n, dt, r, pregunta):
        d = self.calcula_d(n, dt, r, pregunta)
        alpha = np.zeros(self.numptosx)
        beta = np.zeros(self.numptosx)
        Aplus = -1 * r
        Acero = (1 + 2 * r)
        Aminus = -1 * r
        if pregunta == 1:
            alpha[0] = 0
            beta[0] = 1  # viene de n(t, 0) = 1, n(t, 1) = 0
        else:
            alpha[0] = 0
            beta[0] = 0
        for i in range(1, self.numptosx - 1):
            alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
            beta[i] = (d[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)
        return alpha, beta

    def avanza_CN(self, n, nnext, dt, r, pregunta):
        # recibe arreglos n, nnext, alpha y beta
        '''Retorna el valor de la densidad de
        la especie para un paso de tiempo dt, por
        el metodo de Crank-Nicolson'''
        alpha, beta = self.calcula_alpha_y_beta(n, dt, r, pregunta)
        nnext[0], nnext[499] = self.cb
        for i in range(1, self.numptosx - 1):
            nnext[i] = alpha[i] * nnext[i+1] + beta[i]
        return nnext


def condicion_inicial_y_borde(condiciones, numptosx, pregunta, arreglox=None):
    if pregunta == 1:
        ci = np.zeros(len(arreglox))
        for i in range(len(arreglox)):
            ci[i] = np.exp((-arreglox[i]**2)/0.1)
        ci[0], ci[len(arreglox) - 1] = condiciones  # condiciones de borde
    else:
        ci = np.random.uniform(low=-0.3, high=0.3, size=numptosx)
    return ci
'''
# Main Setup P1


numptosx = 500
numptost = 5/dt + 1
xi = np.arange(0, 1, 1/(numptosx - 1))
ci = condicion_inicial_y_borde([1, 0], numptosx, 1, xi)  # n(0,x)
Fisher = Crank_Nicolson([1, 0], ci)
dt = 0.25 * (Fisher.mu/Fisher.gamma) * (Fisher.dx**2)
r = dt*Fisher.gamma/(2*(Fisher.dx**2)*Fisher.mu)

# solucion parte difusion


n = ci
nnext = np.zeros(numptosx)
N = np.zeros((numptosx, int(numptost)))
N[:, 0] = ci
for i in range(1, int(numptost)):
    nnext = Fisher.avanza_CN(n, nnext, dt, r, 1)
    n = nnext
    N[:, i] = nnext

# Graficos


x = np.linspace(0, 1, numptosx)
fig1 = plt.figure(1)
fig1.clf()
for j in range(0, len(N[0]), 100):  # Grafica menos curvas
    plt.plot(x, N[:, j], 'r-')
plt.plot(x, N[:, 0], 'b-')  # Grafica condicion inicial
plt.title(r'Solucion Ecuacion Fisher-KPP ($n(x,t)$) \
para distintos valores de tiempo')
red_patch = mpatches.Patch(color='red', label=r'$n(x,t^*); \ t^* \ fijo $')
blue_patch = mpatches.Patch(color='blue', label=r'$n(x,t=0)$')
plt.legend(handles=[red_patch, blue_patch])
plt.xlabel('Posicion (unidades)')
plt.ylabel(r'Solucion $n(x,t)$ (unidades)')
fig1.savefig('fisherkpp')
plt.grid(True)
plt.show()
'''

# P2 Tarea 6


# Main setup P2


numptosx = 500
np.random.seed(100)
ci = condicion_inicial_y_borde([0, 0], numptosx, 2)
Fisher = Crank_Nicolson([0, 0], ci)
dt = 0.25 * (Fisher.mu/Fisher.gamma) * (Fisher.dx**2)
numptost = 5/dt + 1
r = dt*Fisher.gamma/(2*(Fisher.dx**2)*Fisher.mu)


# solucion parte difusion


n = ci
n[0] = 0
n[499] = 0
nnext = np.zeros(numptosx)
M = np.zeros((numptosx, int(numptost)))
M[:, 0] = ci
for i in range(1, int(numptost)):
    nnext = Fisher.avanza_CN(n, nnext, dt, r, 2)
    n = nnext
    M[:, i] = nnext

# Graficos


x = np.linspace(0, 1, numptosx)
fig2 = plt.figure(2)
fig2.clf()
for j in range(0, len(M[0]), 100):  # Grafica menos curvas
    plt.plot(x, M[:, j], 'r-')
plt.plot(x, M[:, 0], 'b-')  # Grafica condicion inicial
plt.title(r'Solucion Ecuacion Newell-Whitehead-Segel \
($n(x,t)$) para distintos valores de tiempo ($seed=100$)')
red_patch = mpatches.Patch(color='red', label=r'$n(x,t^*); \ t^* \ fijo$ ')
blue_patch = mpatches.Patch(color='blue', label=r'$n(x,t=0)$')
plt.legend(handles=[red_patch, blue_patch])
plt.xlabel('Posicion (unidades)')
plt.ylabel(r'Solucion $n(x,t)$ (unidades)')
fig2.savefig('newell')
plt.grid(True)
plt.show()


# In[ ]:



