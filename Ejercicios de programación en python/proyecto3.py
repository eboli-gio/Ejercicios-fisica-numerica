# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 22:28:38 2020

@author: GIOVANNI
"""
import numpy as np
import math as m
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
def func_init_1(x,L,T_0):   #Definición de funcion en t=0
    if(0<x<L):
        return T_0
    else:
        return 0
def Temperatura(L,t,dx,dt,d,k,c_p,func_init,*args):
    size_x=len(np.arange(0,L+dx,dx))    # Tamaño de valores que guardar en el espacio
    size_t=len(np.arange(0,t+dt,dt))    # Tamaño de valores que guardar en el tiempo
    u=np.zeros((size_t,size_x),dtype=np.float32)    # Matriz de valores que regresar de la temperatura
    alpha=k/(d*c_p)
    step_t=1/(2**16)
    eta=alpha*step_t/(dx**2)
    aux=0.
    for i in range(size_x):
        u[0][i]=func_init(aux, L, T_0)
        aux=aux+dx      
    arr_aux=u[0][:]
    aux_ind=1
    for i in np.arange(step_t,t+step_t,step_t):
         for j in range(1,size_x-1):
             u[aux_ind][j]=arr_aux[j]+eta*(arr_aux[j+1]+arr_aux[j-1]-2*arr_aux[j])
         if(i%dt==0):
            aux_ind=aux_ind+1
            arr_aux=u[aux_ind-1][:]
         else:
            arr_aux=u[aux_ind][:]
    return u 
def Temperatura_anal(dx,dt,n,j,d,c_p):
    alpha=j/(d*c_p)
    u=0
    for k in range(n):
        p=(4*100)/(m.pi*(2*k+1))
        q=(2*k+1)*m.pi/L
        r=(-alpha)*((2*k+1)**2)*(m.pi**2)
        u=u+p*(np.sin(q*(dx))*np.exp(r*dt))
    return u
def grafica(x,y,z,char):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(elev=15,azim=60)
    ax.set_xlabel('t')
    ax.set_ylabel('x')
    ax.set_zlabel('T')
    plt.title(char)
    ax.plot_surface(Y,X,Z,cmap=plt.cm.hot)
    plt.show()
    return 0
def grafica_cont(x,y,z,char):
    plt.title(char)
    plt.contour(x,y,z)
    plt.xlabel('x')
    plt.ylabel('t')
    plt.show()
    return 0
L=1.0
t=0.25
dx=1/32
dt=1/(128)
d=2.699
k=2.37
c_p=0.897
T_0=100
T=Temperatura(L, t, dx, dt, d, k, c_p, func_init_1,T_0)
x=np.arange(0,L+dx,dx)
y=np.arange(0,t+dt,dt)
X,Y=plt.meshgrid(x,y)
Z=Temperatura_anal(X,Y,10,k,d,c_p)
grafica(X,Y,Z,'Gráfica usando la solución analítica')
grafica(X, Y, T,'Gráfica usando la solución numérica')
grafica_cont(X, Y, Z, 'Isotermas usando la solución analítica')
grafica_cont(X, Y, T, 'Isotermas usando la solución numérica')



