# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 22:48:52 2020

@author: GIOVANNI
"""
import numpy as np
import sympy as sp
import matplotlib.pyplot as mplot
def L_i(X,i):
    x=sp.Symbol('x')
    r=0
    l_i=0        
    for k in range(len(X)):
        if k!=i and r==1:
            l_i=l_i*(x-X[k])/(X[i]-X[k])
        if k!=i and r==0:
            l_i=(x-X[k])/(X[i]-X[k])
            r=r+1  
    return l_i
def inter_Lag(X,Y):
    p_Lag=0
    for i in range(len(X)):
        p_Lag=p_Lag+L_i(X,i)*Y[i]
    P=sp.lambdify(x, p_Lag,np)
    return P
def graficar_puntos(X,Y,X_1,Y_1):
    f_1=mplot.plot(X_1,Y_1,'bo')
    f_2=mplot.plot(X,Y,'ro')
    mplot.show(f_1,f_2)
    return 0
def graficar_lineas(X,Y,X_1,Y_1):
    f_1=mplot.plot(X_1,Y_1,'b-')
    f_2=mplot.plot(X,Y,'ro')
    mplot.show(f_1,f_2)
    return 0
X=[0,25,50,75,100,125,150,175,200]
Y=[10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7]
x=sp.Symbol('x')
P_num=inter_Lag(X,Y)
X_1=np.linspace(0, 200,num=40)
Y_1=P_num(X_1)
graficar_lineas(X,Y,X_1,Y_1)
graficar_puntos(X,Y,X_1,Y_1)
X_a=[0,25,50]
Y_a=[10.6,16.0,45.0,]
P_a=inter_Lag(X_a,Y_a)
x_a=np.linspace(0, 50,num=10)
y_a=P_a(x_a)
X_b=[50,75,100]
Y_b=[45.0,83.5,52.8]
P_b=inter_Lag(X_b,Y_b)
x_b=np.linspace(50, 100,num=10)
y_b=P_b(x_b)
X_c=[100,125,150]
Y_c=[52.8,19.9,10.8]
P_c=inter_Lag(X_c,Y_c)
x_c=np.linspace(100,150,num=10)
y_c=P_c(x_c)
X_d=[150,175,200]
Y_d=[10.8,8.25,4.7]
P_d=inter_Lag(X_d,Y_d)
x_d=np.linspace(150,200,num=10)
y_d=P_d(x_d)
Y_t=np.append(y_a,np.append(y_b,np.append(y_c,y_d)))
graficar_lineas(X,Y,X_1,Y_t)
graficar_puntos(X,Y,X_1,Y_t)


