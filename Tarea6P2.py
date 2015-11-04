import numpy as np
import matplotlib.pyplot as plt
#Semilla
seed = 99
np.random.seed(seed)
#Constantes
u = 1.5
y = 0.001
#Condiciones de Borde
nt0 = 0
nt1 = 0
#Datos
Pasos = 500
tTotal = 5         #Tiempo de 0 a 5
xTotal = 1        #x de 0 a 1

dt = tTotal / Pasos   #Paso de tiempo
h = xTotal / Pasos   #Paso de espacio

r = (y * dt) / (2 * h ** 2)   #El cambio de variable
#Vectores vacios para rellenar con datos
n = []              #Junta todos los datos de la solucion
nTotal2 = []     #Junta todas las soluciones y se guarda como '2' para importar con Animaciones
#Crea vectores con ceros para reemplazar con datos
B = np.zeros(Pasos)
a = np.zeros(Pasos)
b = np.zeros(Pasos)
#Aplica las condiciones de borde, incluyendo los numeros al azar
R = np.random.uniform(low=-0.3, high=0.3, size=Pasos)
for i in range(Pasos):
    n.append(R[i])
n[0] = nt0
n[-1] = nt1
#Aplicacion del metodo
N = 1
while N < Pasos:
    for j in range(1, Pasos - 1):                            #Calcula B
        B[j] = (r * n[j+1] + (1-2*r) * n[j] +                #
                r * n[j-1] + n[j] * dt * u * (1 - n[j]*n[j])) #
    Ao = (1 + 2 * r)                                     #Calcula alfa y beta
    A = -1 * r                                             #utilizando los coeficientes
    a[0] = 0                                               #y condiciones de borde
    b[0] = nt0                                            #
    for i in range(1, Pasos):                        #
        a[i] = -A / (Ao + A*a[i-1])                  #
        b[i] = (B[i] - A*b[i-1]) / (A*a[i-1] + Ao)#
    nGuarda = []
    for j in range(len(n)):               #Guarda el actual n 
        nGuarda.append(n[j])         #y calcula el siguiente
    nTotal2.append(nGuarda)       #en reversa
    n[0] = nt0                             #
    n[-1] = nt1                            #
    for i in range(Pasos - 2, 0, -1):#
        n[i] = a[i] * n[i+1] + b[i]      #
    N += 1      #Contador para avanzar con el while
#Plots
x = np.linspace(0, 1, Pasos)
fig = plt.figure()
plt.plot(x, nTotal2[0], label='$t=$ $0.0$')  #Calcula t=0
for i in range(50, 500, 100):                    #Luego t=0.5 para avanzar despues cada 1 seg
        plt.plot(x, nTotal2[i], label='$t=$ $'+str(i*dt)+'$')
plt.legend(loc='lower left')
plt.xlabel('$Posicion$ $en$ $el$ $espacio$ $x$')
plt.ylabel('$Densidad$ $n$')
plt.title('$Densidad$ $v/s$ $posicion,$ $entre$ $t=0.0$ $y$' +
          ' $t=4.5$, $seed=$ $'+str(seed)+'$')
fig.savefig('P2 seed '+str(seed)+'.png')
plt.show()
