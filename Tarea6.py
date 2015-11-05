##############################################################################
'''
Metodos Numericos para la Ciencia e Ingenieria
FI3104-1
Tarea 6
Maximiliano Dirk Vega Aguilera
18.451.231-9
'''
##############################################################################
##############################################################################

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
np.random.seed(184512319)   # fija la semilla para la P2

##############################################################################
##############################################################################
#funciones

def inicializa_n(n, N_pasos, h):
    '''
    Rellena n con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    for i in range(N_pasos):
        x = i * h
        n[i] = np.exp((-x**2)/0.1)
    n[0] = 1
    n[-1] = 0



def calcula_b(b, N_pasos, r, e, mu):
    for j in range(1, N_pasos - 1):
        b[j] = r * n[j+1] + (1-2*r) * n[j] + r * n[j-1] + mu*e*(n[j] - n[j]**2)


def calcula_alpha_y_beta(alhpa, beta, b, r, N_pasos):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 1  # viene de la condicion de borde n(t, 0) = 1
    for i in range(1, N_pasos):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)

def avanza_paso_temporal(n, n_next, alpha, beta, N_pasos):
    n_next[0] = 1
    n_next[-1] = 0
    for i in range(N_pasos - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]

def inicializa_n2(n, N_pasos, h):
    '''
    Rellena n con las condiciones iniciales del problema.
    Se asegura que las condiciones en los bordes sean cero.
    '''
    a = np.random.uniform(low=-0.3, high=0.3, size=N_pasos)
    for i in range(N_pasos):
        n[i] = a[i]
    n[0] = 0
    n[-1] = 0

def calcula_b2(b, N_pasos, r, e, mu):
    for j in range(1, N_pasos - 1):
        b[j] = r * n[j+1] + (1-2*r) * n[j] + r * n[j-1] + mu*e*(n[j] - n[j]**3)


def calcula_alpha_y_beta2(alhpa, beta, b, r, N_pasos):
    Aplus = -1 * r
    Acero = (1 + 2 * r)
    Aminus = -1 * r
    alpha[0] = 0
    beta[0] = 0  # viene de la condicion de borde n(t, 0) = 0
    for i in range(1, N_pasos):
        alpha[i] = -Aplus / (Acero + Aminus*alpha[i-1])
        beta[i] = (b[i] - Aminus*beta[i-1]) / (Aminus*alpha[i-1] + Acero)

def avanza_paso_temporal2(n, n_next, alpha, beta, N_pasos):
    n_next[0] = 0
    n_next[-1] = 0
    for i in range(N_pasos - 2, 0, -1):
        n_next[i] = alpha[i] * n_next[i+1] + beta[i]


##############################################################################
##############################################################################
#Desarrollo P1
'''
dn/dt = g*d2n/dx2 + mn - mn2   n(x,t)
crank-nicolson para el gradiente, euler explicito para mu
0<x<1
gamma = 0.001
mu = 1.5
discretizar el espacio en unos 500 pasos
condiciones:
n(t,0) = 1
n(t,1) = 0
n(0,x) = exp(-x**2 / 0.1)
integrar para almenos t=4
'''
'''
gamma = 0.001
mu = 1.5
N_pasos = 500               #numero de pasos espaciales
h = 1. / (N_pasos-1)        #tamanho del paso espacial #verificar el -1
e = 0.005                   #tamanho del paso espacial
r = gamma * e / (2. * h**2) #cambio de variable para simplificacion del metodo
N_pasos_temporales = int(5. / e)    #numero de pasos en el tiempo


n = np.zeros(N_pasos)
n_next = np.zeros(N_pasos)
b = np.zeros(N_pasos)
alpha = np.zeros(N_pasos)
beta = np.zeros(N_pasos)

inicializa_n(n, N_pasos, h)
n_solucion = np.zeros((N_pasos_temporales, N_pasos))
n_solucion[0, :] = n.copy()

for i in range(1, N_pasos_temporales):
    calcula_b(b, N_pasos, r, e, mu)
    calcula_alpha_y_beta(alpha, beta, b, r, N_pasos)
    avanza_paso_temporal(n, n_next, alpha, beta, N_pasos)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()

x = np.linspace(0, 1, N_pasos)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 20):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(0, 1)

plt.title('Densidad de la especie vs Posicion. A traves del tiempo')
plt.xlabel('Posicion de la especie en una dimension')
plt.ylabel('Densidad de la especie')

plt.show()
plt.draw()
'''

##############################################################################
##############################################################################
#Desarrollo P2
'''
dn/dt = g*d2n/dx2 + mn - mn3   n(x,t)
'''
gamma = 0.001
mu = 1.5
N_pasos = 500               #numero de pasos espaciales
h = 1. / (N_pasos-1)        #tamanho del paso espacial #verificar el -1
e = 0.005                   #tamanho del paso espacial
r = gamma * e / (2. * h**2) #cambio de variable para simplificacion del metodo
N_pasos_temporales = int(5. / e)    #numero de pasos en el tiempo


n = np.zeros(N_pasos)
n_next = np.zeros(N_pasos)
b = np.zeros(N_pasos)
alpha = np.zeros(N_pasos)
beta = np.zeros(N_pasos)

inicializa_n2(n, N_pasos, h)
n_solucion = np.zeros((N_pasos_temporales, N_pasos))
n_solucion[0, :] = n.copy()

for i in range(1, N_pasos_temporales):
    calcula_b2(b, N_pasos, r, e, mu)
    calcula_alpha_y_beta2(alpha, beta, b, r, N_pasos)
    avanza_paso_temporal2(n, n_next, alpha, beta, N_pasos)
    n = n_next.copy()
    n_solucion[i, :] = n.copy()

x = np.linspace(0, 1, N_pasos)

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)

for i in range(0, N_pasos_temporales, 20):
    ax.plot(x, n_solucion[i, :])
ax.set_ylim(-1, 1)

plt.title('Newell-Whitehead-Segel.(seed=184512319)')
plt.xlabel('Posicion en una dimension')
plt.ylabel('n(x,t)')

plt.show()
plt.draw()






##############################################################################
