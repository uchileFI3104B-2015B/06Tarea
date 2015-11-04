#!/usr/bin/env python
# -*- coding: utf-8 -*-

#############################################################################
#                                TAREA 6                                    #
#############################################################################

'''
Universidad de Chile
Facultad de Ciencias Fisicas y Matematicas
Departamento de Fisica
FI3104 Metodos Numericos para la Ciencia y la Ingenieria
Semestre Primavera 2015

Nicolas Troncoso Kurtovic
'''

import matplotlib.pyplot as p
import numpy as np


#############################################################################
#                                                                           #
#############################################################################

'''
Este codigo resuelve una EDP de la forma de la ecuacion del calor mas un 
termino adicional f(m). m es un vecto m[espacio][tiempo]
'''

def reaccion(m, i, n, mu):
  '''
  Devuelve la funcion de reaccion para m_i**n
  '''
  return mu * (m[i][n] - m[i][n]**3)


def rellenar_b(m, i, n, dt, dx, gamma, mu):
  '''
  Retorna el valor de b[i]
  '''
  eta = (gamma * dt) / (2. * dx**2)  

  if i == 0:
    return eta * (m[i+1][n] - m[i][n]) + m[i][n] 
                                       + dt * reaccion(m, i, n, mu)
  
  elif i == (len(m[:][n]) - 1):
    return eta * (-2. * m[i][n] - m[i-1][n]) + m[i][n] 
                                       + dt * reaccion(m, i, n, mu)

  else:
    return eta * (m[i+1][n] - 2 * m[i][n] + m[i-1][n]) + m[i][n] 
                                       + dt * reaccion(m, i, n, mu)


def b_vector(m, n, dt, dx, gamma, mu):
  '''
  Crea un vector del largo de m[:][n] con la forma de la solucion para CN
  combinado con Euler
  '''
  l = len(m[:][n])
  b = np.zeros(l)
  
  for i in range(l):
    b[i] = rellenar_b(m, i, n, dt, dx, gamma, mu)
  
  return b


#############################################################################
#                                                                           #
#############################################################################


m = np.zeros(3,3)
























