# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 13:41:34 2020

@author: GIOVANNI
"""

import scipy.misc as sci_misc
import scipy as aci
import numpy as np
import numpy.linalg as np_al
import sympy as sp
def val_f(f,I):
    v=np.zeros((3,1))
    for i in range(3):
        v[i][0]=f[i][0](I[0][0],I[1][0],I[2][0])
    return v

def val_F(F,I):
    V=np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            V[i][j]=F[i][j](I[0][0],I[1][0],I[2][0])
    return V

def Newton_Raphson(g,G,I):
    X_i=I
    n=0
    while True:
        f=val_f(g,X_i)
        J=val_F(G,X_i)
        J_inv=np_al.inv(J)
        delta_X=-np.matmul(J_inv,f)
        X_i=delta_X+X_i
        n=n+1
        if(delta_X[0][0]<=10**(-6) and delta_X[1][0]<=10**(-6) and delta_X[2][0]<=10**(-6)):
            break
    return X_i


X=np.array([0,25,50,75,100,125,150,175,200])
Y=np.array([10.6,16.0,45.0,83.5,52.8,19.9,10.8,8.25,4.7])
S=np.array([9.34,17.9,41.5,85.5,51.5,21.5,10.8,6.29,4.14])
I=np.array([[63865],[75],[400]])
x,a,b,c,y,s=sp.symbols('x,a,b,c,y,s')#a=a_1=f_r;b=a_2=E_r;c=a_3=gamma^2/4
g=a/((b-x)**2+c)
f_1arg=sp.lambdify([y,s,x], sp.simplify(((y-g)/s**2)*sp.diff(g,a)),np)
f_2arg=sp.lambdify([y,s,x], sp.simplify(((y-g)/s**2)*sp.diff(g,b)),np)
f_3arg=sp.lambdify([y,s,x], sp.simplify(((y-g)/s**2)*sp.diff(g,c)),np)
f_1=np.sum(f_1arg(Y,S,X))
f_2=np.sum(f_2arg(Y,S,X))
f_3=np.sum(f_3arg(Y,S,X))
J_1=sp.lambdify([a,b,c],sp.diff(f_1,a),np)
J_2=sp.lambdify([a,b,c],sp.diff(f_1,b),np)
J_3=sp.lambdify([a,b,c],sp.diff(f_1,c),np)
J_4=sp.lambdify([a,b,c],sp.diff(f_2,a),np)
J_5=sp.lambdify([a,b,c],sp.diff(f_2,b),np)
J_6=sp.lambdify([a,b,c],sp.diff(f_2,c),np)
J_7=sp.lambdify([a,b,c],sp.diff(f_3,a),np)
J_8=sp.lambdify([a,b,c],sp.diff(f_3,b),np)
J_9=sp.lambdify([a,b,c],sp.diff(f_3,c),np)
J=np.array([[J_1,J_2,J_3],[J_4,J_5,J_6],[J_7,J_8,J_9]])
f=np.array([[sp.lambdify([a,b,c],f_1,np)],[sp.lambdify([a,b,c],f_2,np)],[sp.lambdify([a,b,c],f_3,np)]])
sol=Newton_Raphson(f,J,I)
print("Las constante del modelo son: f_r=",sol[0][0],",E_r=",sol[1][0],"y gamma=",2*np.sqrt(sol[2][0]))

