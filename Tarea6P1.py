import numpy as np
import matplotlib.pyplot as plt
#Constantes
u = 1.5
y = 0.001
#Condiciones de Borde
nt0 = 1
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
nTotal1 = []     #Junta todas las soluciones y se guarda como '1' para importar con Animaciones
#Crea vectores con ceros para reemplazar con datos
B = np.zeros(Pasos)
a = np.zeros(Pasos)
b = np.zeros(Pasos)
#Aplica las condiciones de borde, incluyendo la exponencial
for i in range(Pasos):
    x = i * h
    n.append(np.exp(- x ** 2 / 0.1))
n[0] = nt0
n[-1] = nt1
#Aplicacion del metodo
N = 1
while N < Pasos:
    for j in range(1, Pasos - 1):                         #Calcula B
        B[j] = (r * n[j+1] + (1-2*r) * n[j] +             #
                r * n[j-1] + n[j] * dt * u * (1 - n[j]))     #
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
    nTotal1.append(nGuarda)       #en reversa
    n[0] = nt0                             #
    n[-1] = nt1                            #
    for i in range(Pasos - 2, 0, -1):#
        n[i] = a[i] * n[i+1] + b[i]      #
    N += 1      #Contador para avanzar con el while
#Plots
x = np.linspace(0, 1, Pasos)
fig = plt.figure()
for i in range(0, 500, 50):
    if i == 0:      #Destaca el t=0
        plt.plot(x, nTotal1[i], label='         $t=0.0$\n' +
                 '  $Creciendo$ $hacia$\n$la$' +
                 ' $derecha$ $cada$ $0.5$', color='k')
    else:
        plt.plot(x, nTotal1[i], color='b')
plt.legend(loc='lower left')
plt.xlabel('$Posicion$ $en$ $el$ $espacio$ $x$')
plt.ylabel('$Densidad$ $de$ $la$ $especie$ $n$')
plt.title('$Densidad$ $v/s$ $posicion,$ $entre$ $t=0$ $y$ $t=4.5$')
fig.savefig("P1.png")
plt.show()
