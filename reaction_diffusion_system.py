
import numpy as np


class reaction_diffusion_system(object):

    ''' Clase que permite crear ecuaciones de reaccion-difusion
    en 1-D '''

    _n = [] # Almacena vector solucion
    _gamma = 0 # Coeficiente de difusion

    ''' Almacena coeficientes polinomio de reaccion '''
    _reac_poly = [0.0]


    def __init__(self, num_puntos):
        ''' Inicia el vector _n con el tamano solicitado '''

        self._n = np.zeros(num_puntos)

    # END of __init__



    def set_init_conditions(self, vector):
        ''' Setea condiciones iniciales a partir de vector
            Es importante notar que las condiciones iniciales
            deben ser de largo vector -2 para caber considerando
            las Condiciones de borde '''

        len_init_conditions = len(self._n) - 2

        assert len(vector) == len_init_conditions,\
            'Condiciones iniciales de tamano no compatible'

        for i in range(1, len_init_conditions):
            self._n[i] = vector[i]

        print('Condiciones iniciales establecidas')

    # END of set_init_conditions

    def set_border_conditions(self, left, right):
        ''' Escribe las condiciones de borde
        izquierda y derecha del sistema '''

        self._n[0] = left
        self._n[-1] = right

        print('Condiciones de borde fijadas')

    # END of set_border_conditions

    def set_reac_coefficients(self, reac_poly):
        ''' Recibe y almacena los coeficientes del
        polinomio de reaccion de la ecuacion.
        Orden de magnitud va de menor a mayor.

        Ejemplo:
        reac_poly = [2 1 3 -1]
        Componente reaccion: 2 + n + 3n^2 - n^3
        '''

        self._reac_poly = np.array(reac_poly)
        grado_poly = len(self._reac_poly) - 1

        print('Se ha establecido polinomio de reaccion de grado '+ str(grado_poly))

    # END of set_reac_coefficients

    def set_gamma(self, gamma):
        ''' Permite escoger el valor de gamma '''
        self._gamma = gamma

    # END of set_gamma

    def _eval_reac_poly(n):
        ''' Evalua polinomio de reaccion '''
        value = 0
        grade = 0

        for coef in _reac_poly:
            value += coef * n ** grade
            grade += 1

        return value
    # END of _eval_reac_poly

    def _reac_solver(h, n):
        ''' Entrega valor de y_n+1 a partir de h (intervalo temporal)
        e y_n para la componente de reaccion '''
        return n + h * _eval_reac_poly(n)
    # END of _reac_solver

    def integrate(self, start_time, stop_time, step_size):
        ''' Resuelve el sistema desde t = start_time
        hasta stop_time, con resolucion temporal step_size
        Retorna nuevo arreglo con situacion en stop_time '''

        n = np.copy(_n)
        n_next = np.copy(n)
        h = step_size
        t = np.arange(start_time, stop_time + h, h)

        for i in range(t):
            for j in range(1, len(n) - 1):
                n_reac = _reac_solver(h, n[j])

    # END of integrate


def make_reaction_diffusion_system(num_puntos):

    num_puntos = int(num_puntos) # Convierte tamano a entero

    assert num_puntos > 0,\
    'Numero de puntos invalido'

    ''' Permite crear un objeto de sistema reaccion-difusion '''
    new_system = reaction_diffusion_system(num_puntos)

    print('Sistema con ' + str(num_puntos) + ' puntos inicializado')

    return new_system

# END of make_reaction_diffusion_system
