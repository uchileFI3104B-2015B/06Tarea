
'''

'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(456764) #  Fijar semilla


def inicializa_T(T, N_steps, h):
    '''
    Rellena T con las condiciones iniciales del problema.
    '''
    c3 = np.random.uniform(low=-0.3, high=0.3, size=N_steps)

    for i in range(N_steps):
        T[i] = c3[i]

    #  Condiciones iniciales
    T[0] = 0
    T[-1] = 0

#  Crank nicolson y Euler explicito

def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * T[j+1] + (dt * mu * ( 1 - T[j]**2) + 1 - 2 * r) * T[j] + r * T[j-1]


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde T(t, 0) = 1

    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus * alpha[i-1])
        beta[i] = (b[i] - Aminus * beta[i-1]) / (Aminus * alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):
    T_next[0] = 0
    T_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]


# - - - MAIN - - - #

mu = 1.5
gamma = 0.001
t_fin = 4
t_in = 0
N_pasos_temporales = 200
N_steps = 500
numero_pasos_t = N_pasos_temporales - 1
dt = (t_fin - t_in) / (numero_pasos_t)

h = 1 / (N_steps - 1)  # Donde xfin - xini = 1
r = gamma * dt / (2 * h**2)

T = np.zeros(N_steps)
T_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_T(T, N_steps, h)

# Guardar las soluciones en cada paso
T_solucion = np.zeros((N_pasos_temporales, N_steps))
T_solucion[0, :] = T.copy()

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(T, T_next, alpha, beta, N_steps)
    T = T_next.copy()
    T_solucion[i, :] = T.copy()


# - - - PLOTS - - - #

x = np.linspace(0, 1, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 10):
    ax.plot(x, T_solucion[i, :])
ax.set_ylim(-1, 1)
ax.set_xlabel("Cantidad sustancia $x$ [adimensional]")
ax.set_ylabel("Concentracion sustancia $n$ [adimensional]")
ax.set_title("n vs x")


plt.show()
plt.draw()
