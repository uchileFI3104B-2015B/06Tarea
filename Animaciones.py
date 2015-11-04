from Tarea6P1 import nTotal1, x                  #Importa los datos de la pregunta 1 y 2
from Tarea6P2 import nTotal2, seed             #
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

#Crea figura de P1
fig = plt.figure()
Evol = []
for i in range(0, 500, 5):
    Evol.append(plt.plot(x, nTotal1[i], color='b'))
plt.xlabel('$Posicion$ $en$ $el$ $espacio$ $x$')
plt.ylabel('$Densidad$ $de$ $la$ $especie$ $n$')
plt.title('$Densidad$ $v/s$ $posicion,$ $entre$ $t=0.0$ $y$ $t=4.5$')
Evolucion = animation.ArtistAnimation(fig, Evol, interval=25,
                                      repeat_delay=30000, blit=True)
plt.show()

#Crea figura de P2
fig = plt.figure()
Evol = []
for i in range(0, 500, 5):
    Evol.append(plt.plot(x, nTotal2[i], color='b'))
plt.xlabel('$Posicion$ $en$ $el$ $espacio$ $x$')
plt.ylabel('$Densidad$ $de$ $la$ $especie$ $n$')
plt.title('$Densidad$ $v/s$ $posicion, $' +
          ' $entre$ $t=0.0$ $y$ $t=4.5$, $seed=$ $'+str(seed)+'$')
Evolucion = animation.ArtistAnimation(fig, Evol, interval=25,
                                      repeat_delay=30000, blit=True)
plt.show()
