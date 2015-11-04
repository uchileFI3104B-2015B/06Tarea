'''
Este script resuelve la ecuacion de Newell-Whitehead-Segel
(parabolica no lineal) utilizando metodo de Crank-Nicolson para la parte de
difusion y metodo explicito para la parte de reaccion, usando "time stepping".
La ecuacion es de la forma: n_t= gamma n_xx + mu n (1-n^2). Donde n=n(t,x).
Se usa gamma=0.001 y mu=1.5.
Utiliza condiciones de borde u(t,0)=0, u(t,1)=0 y condicion inicial
u(0,x)=np.random.uniform(low=-0.3, high=0.3, size=Nx).
Se resuelve para un tiempo entre 0 y 1, con dt=0.002 y dx=0.002
y se escoge realizar solo 6 plots para obtener una mejor imagen.
Utiliza seed=1876.
'''

from __future__ import division
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.sparse import (spdiags, coo_matrix, csc_matrix,
                          dia_matrix, dok_matrix, identity)


def M(x):
    # me devuelve una matriz de len(x) x len(x), tridiagonal
    # con secuencia 1, -2, 1.
    ones = np.ones(len(x))
    M = spdiags([ones, -2*ones, ones], (-1, 0, 1), len(x), len(x)).tocsc()
    return M


def Imod(x):
    data = np.ones(len(x))
    data[0] = 0
    data[-1] = 0
    offset = 0
    shape = (len(x), len(x))
    Imod = dia_matrix((data, offset), shape=shape).todok()
    return Imod.tocsc()


def solucion_inicial(x, t=0):
    np.random.seed(1876)
    return np.random.uniform(low=-0.3, high=0.3, size=len(x))


# Coeficientes y constantes
x = np.linspace(0, 1, 500)
x_range = (np.max(x) - np.min(x))
dx = x_range / len(x)
gamma = 0.001
mu = 1.5
s = gamma / dx**2
dt = 0.002
valor_izq = 0.0
valor_der = 0.0

# Estado inicial de u
n_inicial = np.matrix(solucion_inicial(x)).reshape((len(x), 1))
n_inicial[0, 0] = valor_izq
n_inicial[-1, -1] = valor_der
n = n_inicial

# Matrices and vectores
M = M(x)                # Matriz de coeficientes
Im = Imod(x)            # matriz identidad modificada
I = identity(len(x))    # matriz identidad
theta = 0.5 * np.ones(len(x))
theta[0] = 0.0
theta[-1] = 0.0
theta_matrix = dia_matrix((theta, 0), shape=(len(x), len(x))).tocsc()
theta_matrix_1m = dia_matrix((1-theta, 0), shape=(len(x), len(x))).tocsc()

STEPS = 2500
PLOTS = 5
fig = plt.figure(1)
fig.clf()
ax = plt.subplot(111)


for i in range(0, STEPS+1):
    # crea una sola lista larga
    n_array = np.asarray(n).flatten()

    # reaccion (parte no lineal, calculada siempre)
    r = np.matrix(mu*n_array*(1.0 - n_array**2)).reshape((len(x), 1))

    if i % abs(STEPS/PLOTS) == 0:
        plot_num = abs(i/(STEPS/PLOTS))
        completado = plot_num / PLOTS
        print "fraccion completada", completado
        print "I: %g" % (simps(n_array, dx=dx), )
        ax.plot(x, n_array, "-o", color=plt.cm.Accent(completado),
                label="t="+str(int(i*dt)))
        ax.legend(loc='center left', bbox_to_anchor=(1., 0.5))

    # Matriz A
    A = (I - dt*s*theta_matrix*M)

    # vector b
    b = csc_matrix((I + dt*s*theta_matrix_1m*M)*n + dt*r)

    # Condiciones de borde
    b[0, 0] = valor_izq
    b[-1, -1] = valor_der

    # Se resuelve An=b
    n = spsolve(A, b)                       # Returns an numpy array,
    n = np.matrix(n).reshape((len(x), 1))   # need a column matrix

plt.subplots_adjust(left=None, bottom=None, right=0.8, top=None,
                    wspace=None, hspace=None)

plt.title("seed=1876, dt=%g." % dt)
plt.xlabel("x")
plt.ylabel("n(x)")
plt.savefig("figura3.png")
plt.show()
plt.draw()
