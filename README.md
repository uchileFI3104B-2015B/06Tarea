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

<img src="img/NWS.png" height="50px"/>

> Latex: `\frac{\partial n}{\partial t} = \gamma \frac{\partial^2n}{\partial x^2} + \mu ( n - n^3)`


Esta vez la ecuación tiene 3 puntos de equilibrio n = 0 (inestable) y n =&pm; 1
(estables). Explique en argumentos simples por qué son estables.

Integre esta ecuación siguiendo la misma estrategia que en la pregunta anterior
(mismas constantes también) pero con las siguientes condiciones de borde:


<img src="img/borde2.png" height="90px"/>

> Latex: `\begin{flalign*} n(t, 0) &= 0\\ n(t, 1) &= 0\\ n(0, x) &= \texttt{np.random.uniform(low=-0.3, high=0.3, size=Nx)} \end{flalign*}`

Si resolvió la pregunta anterior de manera ordenada y modular, entonces sólo
necesitará hacer un par de pequeños cambios a su código.

> Nota: las condiciones iniciales son aleatorias. Asegúrese de setear las
> condiciones de borde (n = 0 para x=0, 1) despues de asignar las condiciones
> aleatorias. También es importante setear la _semilla_ al principio del script
> (`np.random.seed(<algún int>)`), de esa manera su resultado será reproducible
> y no cambiará cada vez que ejecute el script.

Cambie la semilla un par de veces y estudie los cambios en su resultado.

Presente sus resultados mediante los gráficos que le parezcan relevantes e
interprete los resultados.

 __Otras Notas.__

- Utilice `git` durante el desarrollo de la tarea para mantener un historial de
  los cambios realizados. La siguiente [*cheat
  sheet*](https://education.github.com/git-cheat-sheet-education.pdf) le puede
  ser útil. Evaluar el uso efectivo de `git`. Recuerde hacer cambios
  significativos pero relativamente pequeños y guardar seguido.  Evite hacer
  `commits` de código que no compila y deje mensajes que permitan entender los
  cambios realizados.

- Evaluaremos su uso correcto de python. Si define una función relativamente
  larga o con muchos parámetros, recuerde escribir el *doctsring* que describa
  los parametros y que es lo que hace la función.  También recuerde usar nombres
  explicativos para las variables y las funciones.  El mejor nombre es aquel que
  permite entender que hace la función sin tener que leer su implementación.

- Los códigos entregados deben pasar las reglas de
  [PEP8](https://www.python.org/dev/peps/pep-0008/). En la línea de comando
  puede intentar `pep8 <script.py>` y asegurarse de que no hay errores ni
  advertencias. Los comandos `pep8 --show-source <script.py>` y `pep8
  --show-pep8 <script.py>` que vimos en clase le pueden ser útiles. Si es de
  aquellos que resuelven su tarea usando el `ipython notebook`, entonces exporte
  la tarea a un archivo normal de `python` y asegúrese de que pasa el test de
  PEP8.

- La tarea se entrega como un *pull request* en github. El *pull request* debe
  incluir todos los códigos usados además de su informe.

- El informe debe ser entregado en formato *pdf*, este debe ser claro sin
  información ni de más ni de menos.
