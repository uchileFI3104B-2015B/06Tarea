from __future__ import division
import numpy as np
import matplotlib.pyplot as plt


def inicializa_n(n, N_steps, h):
    '''
    Rellena n con las condiciones iniciales t=0 del problema.
    Se consideran las condiciones de borde para x=0 y x=1
    '''
    for i in range(N_steps):
        x = i * h
        n[i] = np.exp(- x ** 2 / 0.1)
    n[0] = borde1
    n[-1] = borde2


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * n[j+1] + (1-2*r) * n[j] + r * n[j-1] + n[j] * (dt * mu - dt * mu * n[j])



def calcula_alpha_y_beta(alhpa, beta, b, r, N_Steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = borde1  # viene de la condicion de borde n(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(n, n_next, alpha, beta, N_steps):
    n_next[0] = borde1
    n_next[-1] = borde2
    for i in range(N_steps - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]

# Main
# setup
mu=1.5
gamma=0.001
N_steps = 500
t_inicial=0
t_final=10
dt=0.01
N_pasos_temporales =((t_final - t_inicial) / dt) + 1
borde1=1
borde2=0
x_inicial=0
x_final=1
h = (x_final-x_inicial) / (N_steps - 1)
r = (gamma*dt) /( 2* h**2)

n = np.zeros(N_steps)
n_next = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)
b = np.zeros(N_steps)


inicializa_n(n, N_steps, h)

# Queremos guardar las soluciones en cada paso
n_solucion = np.zeros((N_pasos_temporales, N_steps))
n_solucion[0, :] = n.copy()

for i in range(1, int(N_pasos_temporales)):
    calcula_b(b, N_steps, r)
    calcula_alpha_y_beta(alpha, beta, b, r, N_steps)
    avanza_paso_temporal(n, n_next, alpha, beta, N_steps)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()


# Grafico p1

x = np.linspace(x_inicial, x_final, N_steps)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, int(N_pasos_temporales), 10):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(0, 1)
ax.set_title('Grafico n v/s x')
ax.set_xlabel('x [unidades arbitrarias]')
ax.set_ylabel('n [unidades arbitrarias]')
plt.savefig('figurap1.png')


plt.show()
plt.draw()
