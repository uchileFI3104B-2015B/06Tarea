'''
Este script...
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def inicializa_T(T, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * h
        T[i] = np.exp(- x ** 2 / 0.1)
    T[0] = 1
    T[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = (T[j+1] * r  + T[j-1] * r +
                T[j] * (1 + e * mu * (1 - T[j])))


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde T(t, 0) = 1
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):
    T_next[0] = 1
    T_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]

# Main

# setup
e = 0.01
gamma = 0.001
mu = 1.5


N_steps = 500
N_pasos_temporales = ( 4 / e ) + 1

h = 1 / (N_steps - 1)

r = gamma * e / (2 * h ** 2)

T = np.zeros(N_steps)
T_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_T(T, N_steps, h)

# Queremos guardar las soluciones en cada paso
T_solucion = np.zeros((N_pasos_temporales, N_steps))
T_solucion[0, :] = T.copy()

for i in range(1, int(N_pasos_temporales)):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(T, T_next, alpha, beta, N_steps)
    T = T_next.copy()
    T_solucion[i, :] = T.copy()

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, int(N_pasos_temporales), 20):
    ax.plot(x, T_solucion[i, :])

ax.set_ylim(0, 1)

plt.show()
plt.draw()
