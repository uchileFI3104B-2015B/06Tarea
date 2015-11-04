'''
Este script...
'''
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(3)

def inicializa_T(T, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    R = np.random.uniform(low=-0.3, high=0.3, size=N_steps)
    for i in range(N_steps):
        x = i * h
        T[i] = R[i]
    T[0] = 0
    T[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = (T[j+1] * r  + T[j-1] * r +
                T[j] * (1 - 2 * r + mu * e * (1 - T[j] * T[j])))


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):
    T_next[0] = 0
    T_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]

# Main

# setup
e = 0.01
gamma = 0.001
mu = 1.5


N_steps = 500
N_pasos_temporales = ( 5 / e ) + 1

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

for i in range(0, int(N_pasos_temporales), 90):
    ax.plot(x, T_solucion[i, :],label="t="+str(i*e))

ax.set_xlabel("Posicion en el espacio x")
ax.set_ylabel("Densidad de la especie n")
ax.set_title("Grafico de densidad versus posicion para distintos tiempos")
plt.legend(loc='lower left')

fig.savefig("graf_2_seed_3.png")
plt.show()
plt.draw()
