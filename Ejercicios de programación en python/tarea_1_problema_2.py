# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 12:55:23 2020

@author: GIOVANNI
"""
import numpy as np 
import tabulate as tab
def esf_bessel_down(x,n,m):
    j=np.zeros(n+m+2,dtype=np.float64)
    j[n+m]=1
    for i in range(n+m+1,0,-1):
        j[i-2]=((2*(i-1)+1)/x)*j[i-1]-j[i]
    j_0=np.sin(x)/x
    res=j[n]*(j_0/j[0])
    return res
def esf_bessel_up(x,n):
    j=np.zeros(n+1,dtype=np.float64)
    j[0]=np.sin(x)/x
    j[1]=np.sin(x)/(x*x)-np.cos(x)/x
    for i in range(n-1):
        j[i+2]=((2*(i+1)+1)/x)*j[i+1]-j[i]
    return j[n]
x=0.1
n=25
l=0
m=50
j_up=np.zeros(3*(n-2)+1,dtype=np.float64)
j_down=np.zeros(3*(n-2)+1,dtype=np.float64)
for i in range(3):
    if l==0:
        for k in range(n-1):
            j_up[k+l]=esf_bessel_up(x,k+3)
            j_down[k+l]=esf_bessel_down(x,k+3,m)
    else:
        for k in range(1,n-1):
            j_up[k+l]=esf_bessel_up(x,k+2)
            j_down[k+l]=esf_bessel_down(x,k+2,m)
    l=l+k-1
    x=x*10
l=0
error=np.zeros(3*(n-2)+1,dtype=np.float64)
for i in range(3):
    if l==0:
        for k in range(n-2):
            error[k+l]=(abs(j_up[k+l]-j_down[k+l]))/(abs(j_up[k+l])+abs(j_down[k+l]))
    else:
        for k in range(1,n-2):
            error[k+l]=(abs(j_up[k+l]-j_down[k+l]))/(abs(j_up[k+l])+abs(j_down[k+l]))
    l=l+k
orden=np.zeros(n-2,dtype=int)
orden[0]=3
for i in range(n-3):
    orden[i+1]=orden[i]+1
tabla1=np.zeros((n-2,4),dtype=np.float64)
for i in range(n-2):
        tabla1[i][0]=orden[i]
        tabla1[i][1]=j_up[i]
        tabla1[i][2]=j_down[i]
        tabla1[i][3]=error[i]
print("Para x=0.1 obtenemos que:")
print(tab.tabulate(tabla1,headers=['l','j_up','j_down','error'],tablefmt='latex',stralign='center'))
tabla2=np.zeros((n-2,4),dtype=np.float64)
for i in range(23,45):
        tabla2[i-23][0]=orden[i-23]
        tabla2[i-23][1]=j_up[i]
        tabla2[i-23][2]=j_down[i]
        tabla2[i-23][3]=error[i]
tabla3=np.zeros((n-2,4),dtype=np.float64)
for i in range(45,68):
        tabla3[i-45][0]=orden[i-45]
        tabla3[i-45][1]=j_up[i]
        tabla3[i-45][2]=j_down[i]
        tabla3[i-45][3]=error[i]
print("Para x=1 obtenemos que:")
print(tab.tabulate(tabla2,headers=['l','j_up','j_down','error'],tablefmt='latex',stralign='center'))
print("Para x=10 obtenemos que:")
print(tab.tabulate(tabla3,headers=['l','j_up','j_down','error'],tablefmt='latex',stralign='center'))


