
import numpy as np
from numpy.linalg import *


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

        print('Gamma inicializado')

    # END of set_gamma

    def _eval_reac_poly(self, n):
        ''' Evalua polinomio de reaccion '''
        value = 0.0
        grade = 0.0

        for coef in self._reac_poly:
            value += coef * (n ** grade)
            grade += 1.0

        return value
    # END of _eval_reac_poly

    def _reac_solver(self, h, n):
        ''' Entrega valor de y_n+1 a partir de h (intervalo temporal)
        e y_n para la componente de reaccion '''
        return h * self._eval_reac_poly(n)
    # END of _reac_solver

    def _diff_solver(self, time_step, x_step, n):
        ''' Integra vector n usando metodo de crank-nicolson '''
        h = time_step
        M = self._gamma * h / (x_step**2)

        # print(len(n))

        len_matrix = len(n) - 2
        A = np.zeros((len_matrix, len_matrix))
        b = np.zeros(len_matrix)

        ''' Se escribe elemento b '''
        for i in range(len_matrix):
            b[i] = M / 2 * n[i] + (1-M) * n[i+1] + M / 2 * n[i+2]

        ''' Bordes de matriz A y vector b '''
        A[0, 0] = (1 + M)
        A[0, 1] = -M / 2
        b[0] = b[0] + M / 2 * n[0]

        for i in range(1, len_matrix - 1):
            A[i, i-1] = -M / 2
            A[i, i] = 1 + M
            A[i, i+1] = -M / 2

        A[-1, -2] = -M / 2
        A[-1, -1] = 1 + M
        b[-1] = b[-1] + M / 2 * n[-1]

        ''' Resolver sist lineal '''
        x = solve(A, b)
        for i in range(0, len(x)):
            n[i+1] = x[i]

        return n
    # END of _diff_solver

    def integrate(self, start_time, stop_time, time_step, x_step):
        ''' Resuelve el sistema desde t = start_time
        hasta stop_time, con resolucion temporal time_step
        Retorna nuevo arreglo con situacion en stop_time
        x_step es el paso en unidad de largo '''

        n = np.copy(self._n)
        n_reac = np.copy(n)
        n_diff = np.copy(n)
        h = x_step
        t = np.arange(start_time, stop_time + h, h)

        for i in range(len(t)):
            n_diff = self._diff_solver(h, x_step, n_diff)

            for j in range(1, len(n) - 1):
                n_reac[j] = self._reac_solver(h, n_reac[j])
                # n[j] = n_diff[j] + n_reac[j]

        for i in range(1,len(n) - 1):
            n[i] = n_diff[i] + n_reac[i]

        return n
    # END of integrate


def make_reaction_diffusion_system(num_puntos):

    num_puntos = int(num_puntos) # Convierte tamano a entero

    assert num_puntos > 0,\
    'Numero de puntos invalido'

    ''' Permite crear un objeto de sistema reaccion-difusion '''
    new_system = reaction_diffusion_system(num_puntos)

    print('Sistema de ' + str(num_puntos) + ' puntos iniciado')

    return new_system

# END of make_reaction_diffusion_system
