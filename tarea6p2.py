'''
Este script resuelve la ecuación de Newell-Whitehead-Segel mediante
el método de Crank-Nicolson.
La ecuación a resolver es:
dn/dt = gamma * d2n/dx2 + mu* (n - n^3)
Condiciones de borde e incial:
n(t, 0) = 0 ; n(t, 1) = 0 ;
n(0, x) = np.random.uniform(low=-0.3, high=0.3, size=Nx)
'''

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

semilla = raw_input("Ingresar semilla:")
np.random.seed(int(semilla))

def inicializa_n(n, N_steps, dx):
    '''
    Rellena n con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    N = np.random.uniform(low=-0.3, high=0.3, size=N_steps)
    for i in range(N_steps):
        n[i] = N[i]
    n[0] = 0
    n[-1] = 0

def calcular_b(b, N_steps, r):
    '''
    Calcula la parte derecha de la ecuación, que depende
    del instante anterior.
    '''
    for j in range(1, N_steps - 1):
        b[j] = r*n[j+1] + ((1 - 2*r) + dt*mu*(1 - n[j]**2))*n[j] + r * n[j-1]
    pass

def alpha_beta(alpha, beta, b, r, N_steps):
    '''
    Calcula los valores de alpha y beta.
    '''
    A_plus = -1 * r
    A_zero = (1 + 2 * r)
    A_minus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde n(t, 0) = 0
    for i in range(1, N_steps):
        alpha[i] = -A_plus / (A_zero + A_minus * alpha[i-1])
        beta[i] = (b[i] - A_minus * beta[i-1]) / (A_minus * alpha[i-1] + A_zero)
    pass

def time_step(n, n_next, alpha, beta):
    '''
    Avanza un instante en el tiempo.
    '''
    n_next[0] = 0 #condición de borde n(t,0) = 0
    n_next[-1] = 0 #condición de borde n(t,1) = 0
    for i in range(N_steps - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]
    pass

x_i = 0
x_f = 1
t_i = 0
t_f = 4
N_steps = 500
N_steps_t = 50
dx = (x_f - x_i) / (N_steps - 1)
dt = (t_f - t_i) / (N_steps_t - 1)
mu = 1.5
gamma = 0.001
r = (gamma * dt) / (2 * dx**2)
n = np.zeros(N_steps)
n_next = np.zeros(N_steps)
b = np.zeros(N_steps)
alpha = np.zeros(N_steps)
beta = np.zeros(N_steps)

inicializa_n(n, N_steps, dx)

n_sol = np.zeros((N_steps_t, N_steps))
n_sol[0, :] = n.copy()

for i in range(1, N_steps_t):
    calcular_b(b, N_steps, r)
    alpha_beta(alpha, beta, b, r, N_steps)
    time_step(n, n_next, alpha, beta)
    n = n_next.copy()
    n_sol[i, :] = n.copy()

fig, ax = plt.subplots()
x = np.linspace(x_i, x_f, N_steps)
ax.set_ylim(0, 1)
ax.set_xlim(0, 1)
ax.set_xlabel("$x$")
ax.set_ylabel("$n$")
for i in range(0, N_steps_t):
    ax.plot(x, n_sol[i, :], color="b")
fig.savefig("grafico2.png")
fig.show()

fig2, ax2 = plt.subplots()
ims = []
ax2.set_xlabel("$x$")
ax2.set_ylabel("$n$")
for add in np.arange(N_steps_t):
    ims.append(plt.plot(x, n_sol[add, :],
               color="b", label="t= " + str(add)))
im_ani = animation.ArtistAnimation(fig2, ims, interval=50,
                                   repeat_delay=3000, blit=True)
plt.show()
