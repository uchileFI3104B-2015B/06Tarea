'''
Este script resuelve numéricamente la ecuación de Fisher-KPP, usando el metodo
de Crank Nicolson y el de Euler explicito. La ecuacion corresponde a
dT/dt = gamma * d2T/dx2 + mu * T - mu * T^3 , con gamma = 0.001 y mu = 1.5 y
condiciones de borde T(t,0) = 0, T(t,1) = 0 y
T(0,x) = np.random.uniform(low=-0.3, high=0.3, size = N_steps)
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(1)


def inicializa_T(T, N_steps):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean CB1 y CB2.
    '''
    R = np.random.uniform(low=-0.3, high=0.3, size=N_steps)
    for i in range(N_steps):
        T[i] = R[i]
    T[0] = CB1
    T[-1] = CB2


# Con solucion Crank Nicolson y Euler explicito
def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = (r * T[j+1] + (1-2*r) * T[j] + r * T[j-1] +
                T[j] * (dt * mu - dt * mu * T[j] * T[j]))


def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = CB1  # viene de la condicion de borde T(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):
    T_next[0] = CB1
    T_next[-1] = CB2
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]


# Main

# Setup
gamma = 0.001
mu = 1.5

# Condiciones de borde para T(t,0) y T(t,1) respectivamente
CB1 = 0
CB2 = 0

N_steps = 500

# Paso temporal dt
dt = 0.01
t_inicial = 0
t_final = 5
N_pasos_temporales = ((t_final - t_inicial) / dt) + 1

x_inicial = 0
x_final = 1
# Paso espacial h
h = (x_final - x_inicial) / (N_steps - 1)


r = (gamma * dt) / (2 * h ** 2)


T = np.zeros(N_steps)
T_next = np.zeros(N_steps)

b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

# Pone las condiciones iniciales
inicializa_T(T, N_steps)

# Queremos guardar las soluciones en cada paso
T_solucion = np.zeros((N_pasos_temporales, N_steps))
T_solucion[0, :] = T.copy()

# Iteracion
for i in range(1, int(N_pasos_temporales)):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(T, T_next, alpha, beta, N_steps)
    T = T_next.copy()
    T_solucion[i, :] = T.copy()

# Plot

x = np.linspace(x_inicial, x_final, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, int(N_pasos_temporales), 80):
    ax.plot(x, T_solucion[i, :],  label="t="+str(i*dt))

ax.set_xlabel("Posicion en el espacio $x$ [adimensional]")
ax.set_ylabel("Densidad $n$ [adimensional]")
ax.set_title("Grafico de densidad versus posicion, entre t=0 y t=4.8")
plt.legend(loc='upper right')

fig.savefig("p_2_seed1.png")
plt.show()
plt.draw()
