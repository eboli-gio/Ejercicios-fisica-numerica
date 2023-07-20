# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 19:50:42 2021

@author: GIOVANNI
"""

from numba import jit
import numpy as np
import math as m
@jit(nopython=True)
def V(j,u):
    V_s=0.0
    d=0.001
    for i in range(0,j):
       V_s=V_s+(u[i+1]**2+u[i]**2)/j
    if(j!=len(u)-1):
        for i in range(j,len(u)-1):
            V_s=V_s+(u[i+1]**2/(i+1)+u[i]**2/i)
    return 2*m.pi*V_s*d
@jit(nopython=True)
def V_ns(u):
    v=np.zeros(len(u))
    for i in range(1,len(v)):
        v[i]=V(i,u)
    return v
@jit(nopython=True)
def eq_10(u,V,E,d):
    for i in range(1,len(u)-1):
        u[i+1]=2*u[i]-u[i-1]+2*d*d*(-2/(i*d)+V[i]-E)*u[i]
    return u
@jit(nopython=True)
def init_1_s(i,d):
    u_1s=np.zeros(i)
    for j in range(i):
        u_1s[j]=4*np.sqrt(2)*j*d*np.exp(-2*j*d)
    return u_1s
@jit(nopython=True)
def init_2_s(i,d):
    u_2s=np.zeros(i)
    for j in range(i):
        r=j*d
        u_2s[j]=1/np.sqrt(2)*(r-(r**2)/2)*np.exp(-r/2)
    return u_2s
@jit(nopython=True)
def init_3_s(i,d):
    u_3s=np.zeros(i)
    for j in range(i):
        r=j*d
        u_3s[j]=2/np.sqrt(27)*(r-(2*r**2)/3+(2*r**3)/27)*np.exp(-r/3)
    return u_3s
@jit(nopython=True)
def init_2_p(i,d):
    u_2p=np.zeros(i)
    for j in range(i):
        r=j*d
        u_2p[j]=1/(2*np.sqrt(6))*(r**2)*np.exp(-r/2)
    return u_2p
@jit(nopython=True)
def init_3_p(i,d):
    u_3p=np.zeros(i)
    for j in range(i):
        r=j*d
        u_3p[j]=(4*np.sqrt(2))/(27*np.sqrt(3))*(r-(r**2)/6)*r*np.exp(-r/3)
    return u_3p
@jit(nopython=True)
def init_3_d(i,d):
    u_3d=np.zeros(i)
    for j in range(i):
        r=j*d
        u_3d[j]=(4)/(81*np.sqrt(30))*(r**3)*np.exp(-r/3)
    return u_3d
def normalizar(u,d):
    N=0.0
    U=np.copy(u)
    k=len(u)-1
    while True:
        if(np.abs(U[k])>np.abs(U[k-1]) and np.abs(U[k-2])>np.abs(U[k-1])):
            break
        else:
            k=k-1
    print(k)
    if k!=1:
        for i in range(len(U)):
            if i>k:
                U[i]=0.0
    for i in range(k-1):
        N=N+(U[i]**2+U[i+1]**2)*(d/2)
    N=np.sqrt(N)
    N=1/N
    U=N*U
    return U
def auto_con(E,E_aux):
    if abs(E-E_aux)<10**-1:
        return False
    else:
        return True
E=2
E_n=2.125
i_max=9000-1
t=0
u_1s=np.zeros(i_max)
u_ns=np.zeros(i_max)
V_1s=np.zeros(i_max)
V_2s=np.zeros(i_max)
dE=-0.1
d=0.001
u_1s=init_1_s(i_max,d)
k=6
if(k==2):
    u_ns=init_2_s(i_max,d)
if(k==3):
    u_ns=init_3_s(i_max,d)
if(k==4):
    u_ns=init_2_p(i_max,d)
if(k==5):
    u_ns=init_3_p(i_max,d)
if(k==6):
    u_ns=init_3_d(i_max,d)
V_1s=V_ns(u_1s)
V_2s=V_ns(u_ns)
u_1s=eq_10(u_1s,V_2s,E,d)
t_01=u_1s[i_max-1]
E_aux=E
E=E+dE
u_1s=normalizar(u_1s,d)
V_1s=V_ns(u_1s)
u_ns=eq_10(u_ns,V_1s,E_n,d)
t_0n=u_ns[i_max-1]
E_naux=E_n
E_n=E_n+dE
u_ns=normalizar(u_ns,d)
V_2s=V_ns(u_ns)
while(auto_con(E,E_aux) and auto_con(E_n,E_naux) or t<50):
      u_1s=eq_10(u_1s,V_2s,E,d)  
      V_1s=V_ns(u_1s)
      t_11=u_1s[i_max-1]
      if t_01*t_11>0:
        E_aux=E
        E=E+dE
        t_01=t_11
      if t_01*t_11<0:
        E_aux=E
        dE=-dE/2
        E=E-dE
        t_01=t_11
      u_ns=eq_10(u_ns,V_1s,E_n,d)
      t_1n=u_ns[i_max-1]
      V_2s=V_ns(u_ns)
      if t_0n*t_1n>0:
        E_naux=E_n
        E_n=E_n+dE
        t_0n=t_1n
      if t_0n*t_1n<0:
        E_naux=E_n
        dE=-dE/2
        E_n=E_n-dE
        t_0n=t_1n
      t=t+1
print(t)
def A(U,d):
    for i in range(1,len(U)-2):
        if i==1:
            I=U[i+2]*U[i+1]+U[i-1]*U[i]+2*U[i]*U[i+1]-2*(U[i+1]**2+U[i]**2)
        else:
            I=U[i+2]*U[i+1]+U[i-1]*U[i]+2*U[i]*U[i+1]-2*(U[i+1]**2+U[i]**2)
    return -1/(4*d)*I
def B(U,d):
    for i in range(1,len(U)-1):
        if i==1:
            I=U[i+1]/(i+1)+U[i]/i
        else:
            I=U[i+1]/(i+1)+U[i]/i
    return -I
A_1=A(u_ns,d)
B_1=B(u_ns,d)
A_n=A(u_1s,d)
B_n=B(u_1s,d)
E_1sns=E_n+A_n+B_n
if(k==2):
   print("El valor de la energía del 1s2s en Hartrees es:",E_1sns,"Y en eV es:",(E_1sns)*13.6044/0.5)
   print("Con un error de entre",(2.146-E_1sns)*100/2.146,"% y",(2.175-E_1sns)*100/2.175,"%")
if(k==3):
    print("El valor de la energía del 1s3s en Hartrees es:",E_1sns,"Y en eV es:",(E_1sns)*13.6044/0.5)
    print("Con un error de entre",(2.061-E_1sns)*100/2.061,"% y",(2.068-E_1sns)*100/2.068,"%")
if(k==4):
    print("El valor de la energía del 1s2p en Hartrees es:",E_1sns,"Y en eV es:",(E_1sns)*13.6044/0.5)
    print("Con un error de entre",(2.124-E_1sns)*100/2.124,"% y",(2.133-E_1sns)*100/2.133,"%")
if(k==5):
    print("El valor de la energía del 1s3p en Hartrees es:",E_1sns,"Y en eV es:",(E_1sns)*13.6044/0.5)
    print("Con un error de entre",(2.055-E_1sns)*100/2.055,"% y",(2.058-E_1sns)*100/2.058,"%")
if(k==6):
    print("El valor de la energía del 1s3p en Hartrees es:",E_1sns,"Y en eV es:",(E_1sns)*13.6044/0.5)
    print("Con un error de entre",(2.0556-E_1sns)*100/2.0556,"% ")