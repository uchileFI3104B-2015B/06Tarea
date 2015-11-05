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
termino adicional f(m). m es un vecto m[espacio,tiempo]
'''

def reaccion(m, i, n, mu):
  '''
  Devuelve la funcion de reaccion para m_i**n
  '''
  return mu * m[i,n] * (1. - m[i,n])


def crear_m_zeros(nx, nt):
  return np.zeros((nx, nt))


def rellenar_b(m, i, n, dt, dx, gamma, mu):
  '''
  Retorna el valor de b[i]
  '''
  eta = (gamma * dt) / (2. * dx**2)  
  if i == 0:
    return eta * (m[i+1,n] - m[i,n]) + m[i,n] + dt * reaccion(m, i, n, mu)
  
  elif i == (len(m[:,0]) - 1):
    return eta * (-2. * m[i,n] + m[i-1,n]) + m[i,n] + dt * reaccion(m, i, n, mu)

  else:
    return eta * (m[i+1,n] - 2 * m[i,n] + m[i-1,n]) + m[i,n] + dt * reaccion(m, i, n, mu)


def b_vector(m, n, dt, dx, gamma, mu):
  '''
  Crea un vector del largo de m[:,n] con la forma de la solucion para CN
  combinado con Euler
  '''
  l = len(m[:,0])
  b = np.zeros(l)
  eta = (gamma * dt) / (2. * dx**2) 
  for i in range(l):
    if i == 0:
      b[i] = 1. + 2. * eta
    else:
      b[i] = rellenar_b(m, i, n, dt, dx, gamma, mu)
  
  return b


def A_matriz(m, gamma, dt, dx):
  '''
  Crea la matriz A del sistema de ecuaciones
  '''
  l = len(m[:,0])
  A = np.zeros((l,l))
  eta = (gamma * dt) / (2. * dx**2) 
  for i in range(l):
    for j in range(l):
      if i == j:
        A[i,j] = 1. + 2. * eta
      if j == (i - 1):
        A[i,j] = - eta
      if j > 1 and j == (i + 1):
        A[i,j] = - eta
  return A      
        

def un_paso(m, n, dt, dx, gamma, mu):
  A = A_matriz(m, gamma, dt, dx)
  b = b_vector(m, n, dt, dx, gamma, mu)
  x = np.linalg.solve(A, b)
  if np.allclose(np.dot(A, x), b):
    for i in range(len(x)):
      try:
        m[i,n+1] = np.copy(x[i])
      except IndexError:
        break
  else:
    print 'Error de resolucion en paso n = ',n


def set_condicion_inicial(m):
  l = len(m[:,0])
  for i in range(l):
    m[i,0] = np.exp( -10. * ( float(i + 1) / float(l) )**2 )
    m[0,i] = 1.


def resolver(nt, nx, gamma, mu):
  dt = 4. / nt
  dx = 1. / nx
  m = crear_m_zeros(nx, nt)
  set_condicion_inicial(m)
  for n in range(len(m[0,:]-1)):
    un_paso(m, n, dt, dx, gamma, mu)
  return m


def graficar1(m):
  '''
  Grafica la matriz V en un diagrama de tonos
  '''
  p.imshow(np.arcsinh(np.transpose(m)) , interpolation = 'nearest')
  extent=([0., 1., 0., 4.])
  p.xlabel('$X$')
  p.ylabel('$T$')
  p.colorbar()
  p.show()   


def graficar2(m):
  '''
  Grafica V en varios subplots
  '''
  fig = p.figure(1)
  fig.clf()
  ax = fig.add_subplot(111)
  lx = len(m[:,0])
  x = np.linspace(0, 1, lx)
  
  for i in range(0, lx, 50):
    ax.plot(x, m[:,i])
  #ax.set_ylim(0,1)
  p.xlabel('$X$')
  p.ylabel('Densidad')
  p.show()

   
#############################################################################
#                                                                           #
#############################################################################

gamma = 0.001
mu = 1.5
nx = 500
nt = 500

m = resolver(nt, nx, gamma, mu)

graficar1(m)

graficar2(m)


