
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
    x_axis = np.linspace(0, 1, num=numero_puntos)
    delta_x = x_axis[1] - x_axis[0]
    t_ini = 0
    t_fin = 4
    delta_t = 0.001

    ''' Crear sistema P2 '''
    sys_p2 = sys.make_reaction_diffusion_system(numero_puntos)

    # Set condiciones de borde
    sys_p2.set_border_conditions(left=0, right=0)

    # Set condiciones iniciales
    np.random.seed(317)
    init_conditions = np.random.uniform(low=-0.3, high=0.3, size=numero_puntos)
    sys_p2.set_init_conditions(init_conditions[1:-1])

    # Set coeficiente de difusion
    sys_p2.set_gamma(0.001)

    # Set componente de reaccion
    mu = 1.5
    reac_coef = [0, mu, 0, -mu]
    sys_p2.set_reac_coefficients(reac_coef)

    # Integrar sistema P2 con semilla 1
    n_1 = sys_p2._n
    n_p2_1 = sys_p2.integrate(t_ini, t_fin, delta_t, delta_x)

    # Integrar para otro valor de semilla
    np.random.seed(119)
    init_conditions = np.random.uniform(low=-0.3, high=0.3,
                                        size=numero_puntos)
    sys_p2.set_init_conditions(init_conditions[1:-1])

    # Integrar sistema P2 con semilla 2
    n_2 = sys_p2._n
    n_p2_2 = sys_p2.integrate(t_ini, t_fin, delta_t, delta_x)

    ''' Graficar '''
    fig = plt.figure(1)
    ax = fig.add_subplot(121)
    plt.plot(x_axis, n_1, label='t=0 [s]')
    plt.plot(x_axis, n_p2_1, 'r', label='t=4 [s]')
    str_title = "Sistema de Newell-Whitehead-Segel\n\
    Tiempo = 4[s], Random density = 317"
    ax.set_title(str_title, y=1.02)
    plt.ylabel('Densidad n', size=14)
    plt.xlabel('Posicion x [m]', size=14)
    plt.grid()
    plt.legend()

    ax = fig.add_subplot(122)
    plt.plot(x_axis, n_2, label='t=0 [s]')
    plt.plot(x_axis, n_p2_2, 'r', label='t=4 [s]')
    str_title = "Sistema de Newell-Whitehead-Segel\n\
    Tiempo = 4[s], Random density = 119"
    ax.set_title(str_title, y=1.02)
    plt.ylabel('Densidad n', size=14)
    plt.xlabel('Posicion x [m]', size=14)
    plt.grid()
    plt.legend(loc=2)

    fig_name1 = "resultados_p2.png"
    plt.savefig(fig_name1)
    plt.tight_layout()
    plt.show()
