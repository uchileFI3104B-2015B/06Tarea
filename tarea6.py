import numpy as np
import matplotlib.pyplot as plt

u = 1.5
y = 0.001

nt0 = 1
nt1 = 0

Pasos = 500
tTotal = 5
xTotal = 1

dt = tTotal / Pasos
h = xTotal / Pasos

r = (y * dt) / (2 * h ** 2)

n = []
nTotal = []

B = np.zeros(Pasos)
a = np.zeros(Pasos)
b = np.zeros(Pasos)

for i in range(Pasos):
    x = i * h
    n.append(np.exp(- x ** 2 / 0.1))
n[0] = nt0
n[-1] = nt1

N=1
while N < Pasos:
    #ActualizaCoeficientes()
    for j in range(1, Pasos - 1):
        B[j] = (r * n[j+1] + (1-2*r) * n[j] +
                r * n[j-1] + n[j] * dt * u * (1 - n[j]))
    Ao = (1 + 2 * r)
    A = -1 * r
    a[0] = 0
    b[0] = nt0  
    for i in range(1, Pasos):
        a[i] = -A / (Ao + A*a[i-1])
        b[i] = (B[i] - A*b[i-1]) / (A*a[i-1] + Ao)
    #avanza_paso_temporal()
    nGuarda = []
    for j in range(len(n)):
        nGuarda.append(n[j])
    nTotal.append(nGuarda)
    n[0] = nt0
    n[-1] = nt1
    for i in range(Pasos - 2, 0, -1):
        n[i] = a[i] * n[i+1] + b[i]
    N += 1
        
x = np.linspace(0, 1, Pasos)

fig = plt.figure()
for i in range(0, 500, 50):
    if i == 0:
        plt.plot(x, nTotal[i], label='         $t=$ $' + str(i*dt) + '$\n  $Creciendo$ $hacia$\n$la$ $derecha$ $cada$ $0.5$', color='k')
    else:
        plt.plot(x, nTotal[i], color='b')
plt.legend(loc='lower left')
plt.xlabel('$Posicion$ $en$ $el$ $espacio$ $x$')
plt.ylabel('$Densidad$ $de$ $la$ $especie$ $n$')
plt.title('$Densidad$ $v/s$ $posicion,$ $entre$ $t=0$ $y$ $t=4.5$')

fig.savefig("P1.png")
plt.show()
plt.draw()
