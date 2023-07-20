# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 12:24:15 2020

@author: GIOVANNI
"""
import scipy.interpolate as sci
import numpy as np
import matplotlib.pyplot as mplot

def graficar_lineas(X,Y,X_1,Y_1):
    f_1=mplot.plot(X_1,Y_1,'b-')
    f_2=mplot.plot(X,Y,'ro')
    mplot.show(f_1,f_2)
    return 0

X=[0,25,50,75,100,125,150,175,200]
Y=[10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7]
P=sci.CubicSpline(X,Y)
x=np.linspace(0, 200,num=40)
y=P(x)
graficar_lineas(X, Y, x, y)
