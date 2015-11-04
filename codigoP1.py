from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def inicializa_T(T, N_steps, dx):
    '''
    Rellena T con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_steps):
        x = i * dx
        T[i] = np.exp((-x**2)/0.1)
    T[0] = 1
    T[-1] = 0


def calcula_b(b, N_steps, r):
    for j in range(1, N_steps - 1):
        b[j] = r * T[j+1] + ((1-2*r) + dt * mu * (1-T[j])) * T[j] + r * T[j-1]


def calcula_alpha_y_beta(alpha, beta, b, r, N_steps):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1
    for i in range(1, N_steps):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)


def avanza_paso_temporal(T, T_next, alpha, beta, N_steps):

    T_next[0] = 1
    T_next[-1] = 0
    for i in range(N_steps - 2, 0, -1):
        T_next[i] = alpha[i] * T_next[i+1] + beta[i]

def plotiar(sol):

    fig, ax = plt.subplots()
    x_values = np.linspace(x_ini, x_fin, largoX)
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Posicion")
    ax.set_ylabel("Densidad")
    ax.set_title("Densidad entre t=0 y t=10")

    for i in range(0, largoT):
        ax.plot(x_values, sol[i, :], color="r")
    fig.savefig("plot1.jpg")

    # esta parte es una animacion que muestra el cambio en la densidad
    fig2, ax2 = plt.subplots()
    ims = []
    ax2.set_xlabel("Posicion")
    ax2.set_ylabel("Densidad")
    ax2.set_title("Densidad entre t=0 y t=10")
    for add in np.arange(largoT):
        ims.append(plt.plot(x_values, sol[add, :], color="b", label="t= " + str(add)))
    im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
                                       blit=True)
    plt.show()


#Implementamos le codigo


t_ini = 0
t_fin = 10
x_ini = 0
x_fin = 1

largoT = 100
largoX = 500

dt = (t_fin - t_ini) / (largoT - 1)
dx = (x_fin - x_ini) / (largoX - 1)

gamma = 0.001
r = (gamma * dt) / (2 * dx**2)

mu = 1.5

T = np.zeros(largoX)
T_next = np.zeros(largoX)

b = np.zeros(largoX)
alpha = np.zeros(largoX)
beta = np.zeros(largoX)

inicializa_T(T, largoX, dx)


# Queremos guardar las soluciones en cada paso
T_solucion = np.zeros((largoT, largoX))
T_solucion[0, :] = T.copy()

for i in range(1, largoT):
    calcula_b(b, largoX, r)
    calcula_alpha_y_beta(alpha, beta, b, r, largoX)
    avanza_paso_temporal(T, T_next, alpha, beta, largoX)
    T = T_next.copy()
    T_solucion[i, :] = T.copy()

plotiar(T_solucion)
