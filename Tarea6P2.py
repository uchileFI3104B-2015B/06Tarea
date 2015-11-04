import numpy as np
import matplotlib.pyplot as plt

seed = 99
np.random.seed(seed)

u = 1.5
y = 0.001

nt0 = 0
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

R = np.random.uniform(low=-0.3, high=0.3, size=Pasos)

for i in range(Pasos):
    n.append(R[i])
n[0] = nt0
n[-1] = nt1

N=1
while N < Pasos:
    for j in range(1, Pasos - 1):
        B[j] = (r * n[j+1] + (1-2*r) * n[j] +
                r * n[j-1] + n[j] * dt * u * (1 - n[j]*n[j]))
    Ao = (1 + 2 * r)
    A = -1 * r
    a[0] = 0
    b[0] = nt0  
    for i in range(1, Pasos):
        a[i] = -A / (Ao + A*a[i-1])
        b[i] = (B[i] - A*b[i-1]) / (A*a[i-1] + Ao)
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
plt.plot(x, nTotal[0], label='$t=$ $0.0$') 
for i in range(50, 500, 100):
          plt.plot(x, nTotal[i], label='$t=$ $'+str(i*dt)+'$')    
plt.legend(loc='lower left')
plt.xlabel('$Posicion$ $en$ $el$ $espacio$ $x$')
plt.ylabel('$Densidad$ $de$ $la$ $especie$ $n$')
plt.title('$Densidad$ $v/s$ $posicion,$ $entre$ $t=0.0$ $y$ $t=4.5$, $seed=$ $'+str(seed)+'$')
fig.savefig('P2 seed '+str(seed)+'.png')
plt.show()
