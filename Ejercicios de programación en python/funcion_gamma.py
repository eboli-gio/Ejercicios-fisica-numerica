# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 15:56:36 2020

@author: GIOVANNI
"""
import sympy as sp
from scipy.special import gamma
import numpy as np
def gamma_trapezoide(x):
    a=0
    b=2**9
    sum_f=0.0
    steps=np.int64(2**16)
    u=np.linspace(a,b,num=steps)
    h=(b-a)/len(u)
    t=sp.symbols('t')
    f=t**(x-1)*sp.exp(-t)
    gamma=sp.lambdify(t, f,"numpy")
    sum_f=gamma(u[0])+gamma(u[len(u)-1])
    for i in range(1,len(u)-1):
        sum_f=sum_f+2*gamma(u[i])
    sol=(h/2)*sum_f
    return sol

def gamma_simpson(x):
    a=0
    b=2**9
    sum_f=0.0
    steps=np.int64(2**16)
    u=np.linspace(a,b,num=steps)
    t=sp.symbols('t')
    f=t**(x-1)*sp.exp(-t)
    h=(b-a)/steps
    gamma=sp.lambdify(t, f,"numpy")
    for i in range(0,len(u)-2,2):
       sum_f=sum_f+(h/3)*(gamma(u[i])+4*gamma(u[i+1])+gamma(u[i+2]))
           
    return sum_f


x=sp.Rational(3,2)    

print("El valor de la función gamma en 3/2 en scipy es:",gamma(1.5))
print("El valor de la función gamma en 3/2 por trapecio es:",gamma_trapezoide(x))
print("El valor de la función gamma en 3/2 por Simpson es:",gamma_simpson(x))


