
import numpy as np
import matplotlib.pyplot as plt

'''
uso: main.py
'''

import reaction_diffusion_system as sys

# Main driver
if __name__ == '__main__':

    np.random.seed(317)

    numero_puntos = 500.0
    largo = 1.0
    x_axis = np.linspace(0, 1, num = numero_puntos)
    delta_x = x_axis[1] - x_axis[0]
    t_ini = 0
    t_fin = 4
    delta_t = 0.001

    ''' Crear sistema P2 '''
    sys_p2 = sys.make_reaction_diffusion_system(numero_puntos)

    # Set condiciones de borde
    sys_p2.set_border_conditions(left = 0, right = 0)

    # Set condiciones iniciales
    init_conditions = np.random.uniform(low = -0.3, high = 0.3, size = numero_puntos)
    sys_p2.set_init_conditions(init_conditions[1:-1])

    # Set coeficiente de difusion
    sys_p2.set_gamma(0.001)

    # Set componente de reaccion
    mu = 1.5
    reac_coef = [0, mu, 0, -mu]
    sys_p2.set_reac_coefficients(reac_coef)

    # Integrar sistema P2
    n_p2 = sys_p2.integrate(t_ini, t_fin, delta_t, delta_x)

    plt.plot(x_axis, sys_p2._n,'b')
    plt.plot(x_axis, n_p2, 'r')
    plt.show()
