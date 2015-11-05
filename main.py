
import numpy as np

'''
uso: main.py
'''

import reaction_diffusion_system as sys

# Main driver
if __name__ == '__main__':

    numero_puntos = 500.0
    largo = 1.0
    x_axis = np.linspace(0, 1, num = numero_puntos)


    ''' Crear sistema P1 '''
    sys_p1 = sys.make_reaction_diffusion_system(numero_puntos)

    # Set condiciones de borde
    sys_p1.set_border_conditions(left = 1, right = 0)

    # Set condiciones iniciales
    init_conditions = np.exp(-x_axis * x_axis / 0.1)
    sys_p1.set_init_conditions(init_conditions[1:-1])

    # Set coeficiente de difusion
    sys_p1.set_gamma(0.001)

    # Set componente de reaccion
    mu = 1.5
    reac_coef = [0, mu, -mu]
    sys_p1.set_reac_coefficients(reac_coef)
