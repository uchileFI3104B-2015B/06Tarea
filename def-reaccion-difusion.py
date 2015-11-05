
import numpy as np


class box(object):

    ''' Clase que permite crear ecuaciones de reaccion-difusion
    en 1-D '''

    _n = [] # Almacena vector solucion
    _gamma = 0 # Coeficiente de difusion

    ''' Almacena coeficientes polinomio de reaccion '''
    _reac_poly = []


    def __init__(self, num_puntos):
        ''' Inicia el vector _n con el tamano solicitado '''

        self._n = np.zeros(num_puntos)

    # END of __init__


    def set_init_conditions(self, vector):
        ''' Setea condiciones iniciales a partir de vector
            Es importante notar que las condiciones iniciales
            deben ser de largo vector -2 para caber considerando
            las Condiciones de borde '''

        assert len(vector) == (len(self._n) - 2),\
            'Condiciones iniciales de tamano no compatible'

        for i in range(1, len(self._n)):
            self._n[i] = vector[i]

    # END of set_init_conditions

    def set_border_conditions(self, left, right):
        ''' Escribe las condiciones de borde
        izquierda y derecha del sistema '''

        self._n[0] = left
        self._n[-1] = right

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

    # END of set_reac_coefficients

    def set_gamma(self, gamma):
        ''' Permite escoger el valor de gamma '''
        self._gamma = gamma

    # END of set_gamma
