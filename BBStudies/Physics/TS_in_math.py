#!/usr/bin/python3

import numpy as np
import scipy.integrate as integrate

def g(t,r):
    return np.sqrt(1+(r**2-1)*t)
#---------------------------------------
def Q0z(t,azbar,dzbar):
    X =  t*azbar*dzbar
    Y = -t*azbar**2/4
    return np.exp(-t/2*(azbar-dzbar)**2)*Bess2D(X,Y,0)
#---------------------------------------
#---------------------------------------
def dQ0daz(t,azbar,dzbar,etaz):
    X    =  t*azbar*dzbar
    Y    = -t*azbar**2/4
    return etaz*np.exp(-t/2*(azbar-dzbar)**2)*(-azbar/2*(Bess2D(X,Y,0)+Bess2D(X,Y,2)) + dzbar*Bess2D(X,Y,1))
#---------------------------------------

def Bess2D_generating_real(phi,X,Y,n):
    arg =  - X*np.sin(phi) + 2*Y*np.sin(phi)**2
    return np.cos(n*phi)*np.exp(arg)
 
def Bess2D_generating_imag(phi,X,Y,n):
    arg =  - X*np.sin(phi) + 2*Y*np.sin(phi)**2
    return -np.sin(n*phi)*np.exp(arg)

def Bess2D(X,Y,n):
 
    generatingFun,sign = {0:(Bess2D_generating_real,(-1)**(n/2)),
    1:(Bess2D_generating_imag,(-1)**((n-1)/2))}[n%2]
     
    integratedFactor = integrate.quad(lambda phi: generatingFun(phi,X,Y,n), 0, 2*np.pi)[0]

    return sign*np.exp(-X-2*Y)/2/np.pi * integratedFactor

#---------------------------------------
def dC00dx_generating(t,r,ax,ay,dx,dy):
    # bar variables
    axbar,aybar,dxbar,dybar = ax*r, ay/g(t,r) , dx , r*dy/g(t,r)

    # modified version of eta to include all the prefactors
    etax_modif = r/g(t,r)

    return dQ0daz(t,axbar,dxbar,etax_modif)*Q0z(t,aybar,dybar)

def dC00dy_generating(t,r,ax,ay,dx,dy):
    # bar variables
    axbar,aybar,dxbar,dybar = ax*r, ay/g(t,r) , dx , r*dy/g(t,r)

    # modified version of eta to include all the prefactors
    etay_modif = 1/(g(t,r)**2)

    return dQ0daz(t,aybar,dybar,etay_modif)*Q0z(t,axbar,dxbar)


def dC00dx(r,ax,ay,dx,dy):
    return integrate.quad(lambda t: dC00dx_generating(t,r,ax,ay,dx,dy), 0, 1)[0]

def dC00dy (r,ax,ay,dx,dy):
    return integrate.quad(lambda t: dC00dy_generating(t,r,ax,ay,dx,dy), 0, 1)[0]


m = 4
k = 0


y = np.loadtxt("par.in")
array = np.array(y)
dx, dy , r, ax , ay =array
print(dx, dy , r,    ax ,  ay)

print(2*dC00dx(r,ax,ay,dx,dy)/ax)
print(2*dC00dy(r,ax,ay,dx,dy)/ay)