# Tarea nro. 6
## FI3104B - Métodos Numéricos para la Ciencia y la Ingeniería
#### Prof. Valentino González

#P1

La ecuación de Fisher-KPP es una llamada ecuación de reacción-difusión que busca
modelar el comportamiento de una especie animal.  A continuación se presenta su
versión en 1D:

<img src="img/fisher-kpp.png" height="50px"/>

> Latex: `\frac{\partial n}{\partial t} = \gamma \frac{\partial^2n}{\partial x^2} + \mu n - \mu n^2`


La variable n = n(t, x) describre la densidad de la especie como función del
tiempo y la posición. Los 3 términos del lado derecho corresponden a:

- &mu;n : la tendencia de la especie a crecer indefinidamente (suponiendo que
  los recursos sean infinitos).
- –&mu;n<sup>2</sup> : Despues de un tiempo, el aumento en densidad creará
  competencia por los recursos, lo que tenderá a disminuir la densidad.
- &gamma; &nabla; n : La tendencia de la especie a dispersarse para encontrar
  más recursos.
  
La ecuación tiene dos puntos de equilibrio n=0 y n=1, pero solo el segundo es
estable. Las soluciones tienen un comportamiento que es una mezcla de difusión y
un pulso viajero.

Para resolver la ecuación discretice la parte de difusión usando el método de
Crank–Nicolson, y el método de Euler explícito para la parte de reacción.
Resuelva la ecuación para x entre 0 y 1 con &gamma; = 0.001 y &mu; = 1.5.
Discretice el espacio en aproximadamente 500 puntos y considere las siguientes
condiciones de borde:

<img src="img/borde1.png" height="90px"/>

> Latex: `\begin{flalign*} n(t, 0) &= 1\\ n(t, 1) &= 0\\ n(0, x) &= e^{-x^2/0.1} \end{flalign*}`

Por último, elija su paso temporal de modo que la solución sea estable e integre
hasta al menos t = 4 (en las unidades en que están escritas las ecuaciónes y las
constantes).

Presente la solución encontrada e interprete los resultados.


#P2

La ecuación de Newell-Whitehead-Segel es otra ecuación de reacción-difusión que
describe fenómenos de convección y combustión entre otros. La ecuación es la
siguiente:


