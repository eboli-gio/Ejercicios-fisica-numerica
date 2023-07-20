# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 19:45:44 2020

@author: GIOVANNI
"""

import sympy as sp
from scipy.special import ellipk ,ellipe  
import numpy as np

def trapecio_elip(a,b,x,f):
    h=10**(-5)
    n=np.int64(np.round((b-a)/h,decimals=0))
    P=np.linspace(a,b ,num=n)
    sol=f(x,P[0])+f(x,P[len(P)-1])
    for i in range(1,len(P)-1):
        sol=sol+2*f(x,P[i])
    sol=(h/2)*sol
    return sol
def simpson_elip(a,b,x,f):
    h=10**(-5)
    n=np.int64(np.round((b-a)/h,decimals=0))
    P=np.linspace(a,b ,num=n)
    sum_f=0.0
    for i in range(0,len(P)-2,2):
       sum_f=sum_f+(h/3)*(f(x,P[i])+4*f(x,P[i+1])+f(x,P[i+2])) 
    return sum_f

x,t=sp.symbols('x,t')
F=1/sp.sqrt(1-(x*sp.sin(t))**2,[x,t])
F_num=sp.lambdify([x,t], F, np)
G=sp.sqrt(1-(x*sp.sin(t))**2,[x,t])
G_num=sp.lambdify([x,t],G,np)
K_t=trapecio_elip(0,np.pi/2,0.5,F_num)
K_s=simpson_elip(0,np.pi/2,0.5,F_num)
K=ellipk(0.25)
print("La integral elíptica de primera especie valuado en 1/2 por trapecio es:",K_t)
print("La integral elíptica de primera especie valuado en 1/2 por Simpson es",K_s)
print("La integral elíptica de primera especie valuado en 1/2 por scipy es",K)
E_t=trapecio_elip(0,np.pi/2,0.5,G_num)
E_s=simpson_elip(0,np.pi/2,0.5,G_num)
E=ellipe(0.25)
print("La integral elíptica de segunda especie valuado en 1/2 por trapecio es:",E_t)
print("La integral elíptica de segunda especie valuado en 1/2 por Simpson es:",E_s)
print("La integral elíptica de segunda especie valuado en 1/2 por scipy es:",E)

