#VERSAO JAN 2018

import numpy as np
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import *

#equacao de movimento para x2pontos
def Fx(m,x,y):
    return (-x-2*x*y)/m

#equacao de movimento para y2pontos
def Fy(m,x,y):
    return (-x**2 -y +y**2)/m

#mesmas equacoes em coordenadas polares
def Fr(m,r,th): #r2pontos
    return (-r -3*(r**2)*np.sin(th)*(np.cos(th))**2 - (r**2)*np.sin(th)**3)/m

def Fth(m,r,th): #theta2pontos
    return (r**3)*(2*np.cos(th)*(np.sin(th))**2 - np.cos(th)**3 + np.cos(th)*np.sin(th)**2)/m

#Hamiltoniana
def H(vx,vy,x,y):
    return 0.5*(vx**2 + vy**2 + x**2 + y**2) + (x**2)*y -(y**3)/3

#Potencial
def V(x,y):
    return 0.5*(x**2 + y**2) + y*(x**2) - (y**3)/3

#Hamiltoniana em Coordenadas Polares
def Hp(r,theta,vr,vth):
    return 0.5*(vr**2 + r**2) + r*vr*(np.sin(theta)*np.sin(vth) + np.cos(theta)*np.cos(vth)) + (1/2)*r**2 + (np.sin(theta)*np.cos(theta)**2 - (1/3)*np.sin(theta)**3)*r**3
#return 0.5*(vr**2 + vth**2) + 0.5*(r**2) + (r**3)*(np.sin(theta)- 4/3*np.sin**3(theta))

def Pot(a,b):
    return (1/2)*a**2+(np.sin(b)*(np.cos(b)**2)-(1/3)*np.sin(b)**3)*a**3

# Velocidade angular, obtem-se de Hp fazendo vr = 0
"""
def vth(c,r,theta):
    if (r != 0):
        return  2*c/r -2*r*(np.sin(theta)-4/3*(np.sin(theta))**3) -1
"""

# H para K = 0, formando surperficies equipotenciais
N = 5000
h = 5/N #tamanho do passo
yl = np.linspace(-1.5,1.5,N) #
z = list(range(N))
c = [i* .4/20 for i in range(1,20)] #conjunto de valores para energia potencial

def xl(yl,c):
    return (c + 1/3*yl**3 - 1/2*yl**2)/(1/2 + yl)

fig = plt.figure()
for hl in c: #
    if (hl <= 0.16):
        plt.plot(np.sqrt(xl(yl,hl)),yl,'r-')
        plt.plot(-np.sqrt(xl(yl,hl)),yl,'r-')
    else:
        plt.plot(np.sqrt(xl(yl,hl)),yl,'k-')
        plt.plot(-np.sqrt(xl(yl,hl)),yl,'k-')
plt.axis([-1.0,1.0,-1.0,1.0])
plt.show()

#definicao de posicoes x, y e r e theta
x = list(range(N))
y = list(range(N))
r = list(range(N))
th = list(range(N))

#definicao de velocidades vx, vy e vr e vtheta
vx = list(range(N))
vy = list(range(N))
vr = list(range(N))
vth = list(range(N))

#definicao dos termos do Runge-Kutta de 4a ordem (em cartesiano)
k1vx = list(range(N))
k1vy = list(range(N))
k1x = list(range(N))
k1y = list(range(N))

k2vx = list(range(N))
k2vy = list(range(N))
k2x = list(range(N))
k2y = list(range(N))

k3vx = list(range(N))
k3vy = list(range(N))
k3x = list(range(N))
k3y = list(range(N))

k4vx = list(range(N))
k4vy = list(range(N))
k4x = list(range(N))
k4y = list(range(N))

#definicao dos termos do Runge-Kutta de 4a ordem (em polares)
k1vr = list(range(N))
k1vth = list(range(N))
k1r = list(range(N))
k1th = list(range(N))

k2vr = list(range(N))
k2vth = list(range(N))
k2r = list(range(N))
k2th = list(range(N))

k3vr = list(range(N))
k3vth = list(range(N))
k3r = list(range(N))
k3th = list(range(N))

k4vr = list(range(N))
k4vth = list(range(N))
k4r = list(range(N))
k4th = list(range(N))

#definicao das condicoes iniciais de r e theta
a=np.linspace(0,0.6,10) #para aplicar ao r
b=np.linspace(0,2*pi,10) #para aplicar ao theta
rCond=[]
thCond=[]

for i in range(0,9):
    for j in range(0,9):
        if ( Pot(a[i],b[j])<0.16 and Pot(a[i],b[j])>0.0):
            print(Pot(a[i],b[j]))
            rCond.append(a[i])
            thCond.append(b[j])

#print('rCond:',rCond)
#print('thCond:',thCond)
#print('a:',a)
#print('b:',b)
#print(len(rCond),len(a))
#print(len(thCond),len(b))

#Fazendo Runge kutta nas equações de movimento para r e theta ARRUMAR COND INICIAIS A PARTIR DE H !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
d = np.linspace(0,pi,100) #condicoes iniciais de theta
f = np.linspace(0,0.4,100) #condicoes iniciais de r

for j in range(0,len(rCond)): #for j in range(0,len(d)):
    r[0] = rCond[j] #r[0] = f[j]
    vr[0] = 0#velocidade radial
    for k in range(0,len(thCond)):
        th[0] = thCond[k]
        aux = 2*0.17 - r[0]**2 - 2*(r[0]**3)*(np.sin(th[0]) - 4/3*(np.sin(th[0]))**3)
        vth[0] = np.sqrt(aux) #velocidade angular vth[0] = np.sqrt(aux)
        for i in range(0,N-1): #for i in range(0,N-1):
            k1vr[i] = Fr(1,r[i],th[i])*h
            k1vth[i] = Fth(1,r[i],th[i])*h
            k1r[i] = vr[i]*h
            k1th[i] = vth[i]*h

            k2vr[i] = Fr(1,r[i]+k1r[i]/2,th[i]+k1th[i]/2)*h
            k2vth[i] = Fth(1,r[i]+k1r[i]/2,th[i]+k1th[i]/2)*h
            k2r[i] = (vr[i] + k1vr[i]/2)*h
            k2th[i] = (vth[i] + k1vth[i]/2)*h

            k3vr[i] = Fr(1,r[i]+k2r[i]/2,th[i]+k2th[i]/2)*h
            k3vth[i] = Fth(1,r[i]+k2r[i]/2,th[i]+k2th[i]/2)*h
            k3r[i] = (vr[i] + k2vr[i]/2)*h
            k3th[i] = (vth[i] + k2vth[i]/2)*h

            k4vr[i] = Fr(1,r[i]+k3r[i],th[i]+k3th[i])*h
            k4vth[i] = Fth(1,r[i]+k3r[i],th[i]+k3th[i])*h
            k4r[i] = (vr[i] + k3vr[i])*h
            k4th[i] = (vth[i] + k3vth[i])*h

            vr[i+1] = vr[i] + 1/6*(k1vr[i] + 2*(k2vr[i] + k3vr[i]) + k4vr[i])
            vth[i+1] = vth[i] + 1/6*(k1vth[i] + 2*(k2vth[i] + k3vth[i]) + k4vth[i])

            r[i+1] = r[i] + 1/6*(k1r[i] + 2*(k2r[i] + k3r[i]) + k4r[i])
            th[i+1] = th[i] + 1/6*(k1th[i] + 2*(k2th[i] + k3th[i]) + k4th[i])

            x[i+1]=r[i+1]*np.cos(th[i+1])
            y[i+1]=r[i+1]*np.sin(th[i+1])

#    print(x[0],x[N-1])
#    print(r)
#    print(np.cos(th))
#    print(len(x))
#    print(type(r))
#    print(type(th))
#        if (x[i]>pow(10,-1) or y[i]>pow(10,-1)):
#            print(x[i],y[i])
#    if(x[100]> 0 and y[100] > 1):
#        plt.plot(x,y,'g-')
#    if(x[100]> sqrt(3)/2 and y[100] > 0.5):
#        plt.plot(x,y,'b-')
#    if(x[100] < sqrt(3)/2 and y[100] < 0.5):
#        plt.plot(x,y,'y-')
#    else:
#        plt.plot(x,y,'p-')
    plt.plot(x,y)

plt.axis([-1.5,1.5,-1.5,1.5])
plt.plot(np.sqrt(xl(yl,0.16)),yl,'r-')
plt.plot(-np.sqrt(xl(yl,0.16)),yl,'r-')
plt.show()
