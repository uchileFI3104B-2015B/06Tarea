
import numpy as np
import matplotlib.pyplot as plt

'''
uso: main.py
'''

import reaction_diffusion_system as sys

# Main driver
if __name__ == '__main__':

    numero_puntos = 500.0
    largo = 1.0
    x_axis = np.linspace(0, 1, num = numero_puntos)
    delta_x = x_axis[1] - x_axis[0]
    t_ini = 0
    t_fin = 4
    delta_t = 0.001

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
    # reac_coef = [0]
    sys_p1.set_reac_coefficients(reac_coef)

    n_p1 = sys_p1.integrate(t_ini, t_fin, delta_t, delta_x)

    plt.plot(x_axis, sys_p1._n,'r')
    plt.plot(x_axis, n_p1)
    plt.show()
