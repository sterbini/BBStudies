from numpy import *
import scipy.integrate as integrate
import numpy as np

"""
Compute Fourier Coefficient Dmk, and   tune-shifts DQX,Y
for a long-range collision (BBLR),  or an ideal wire.
All functions expect a single value for all parameters.
"""
gflag=0
def set_global(g):
         global gflag
         gflag=g

def g(t,r):
         return sqrt(1+(r**2-1)*t)

# 2D-Bessel LAMDA function

def L2(U1,U2,n):
    return exp(-U1-U2)*integrate.quad(lambda phi:
                cos(n*phi) * exp( U1*cos(phi) + U2*cos(2*phi)),
                -pi/2, 3/2*pi)[0]/2/pi

# Bar-kernel nominator x, or y component Q_z (z=x,y)
# bar theta = (axb, dxb, ayb, dyb)

def QmzL(m,azb,dzb):
    U1 =  azb*dzb
    U2 = -azb**2/4
#    return exp(-1/2*(azb-dzb)**2) * L2(U1,U2,m)
    if abs(azb-dzb)<10:
            return exp(-1/2*(azb-dzb)**2) * L2(U1,U2,m)
    else:
             return 0.0

# BBLR bar-kernel nominator 

def QmkL(m,k,axb,ayb,dxb,dyb):   
    return QmzL(m,axb,dxb)*QmzL(k,ayb,dyb)

# BBLR real-valued Fourier Coefficient 
def Dmk(m,k,ax,ay,dx,dy,r):   
    psix=dx/ax/r
    psiy=ay/ax/r
    psiz=dy/ax
    _g   = lambda xi: g((xi/r/ax)**2,r)
    _i   = 0
    _f   = r*ax
    return 2*integrate.quad(lambda xi: 
    1/xi/_g(xi)*QmkL(m=m, k=k, 
                axb=xi,
                ayb=xi*psiy/_g(xi),
                dxb=xi*psix,
                dyb=xi*psiz/_g(xi)), _i, _f)[0]

# Tune-shift bar-kernel nominator 

def TXL( axb,ayb, dxb,dyb):
    return QmzL(0,ayb,dyb)*(
               QmzL(0,axb,dxb)+QmzL(2,axb,dxb)-2*dxb/axb* QmzL(1,axb,dxb)
                                         )
def TYL( axb,ayb, dxb,dyb):
    return QmzL(0,axb,dxb)*(
               QmzL(0,ayb,dyb)+QmzL(2,ayb,dyb)-2*dyb/ayb* QmzL(1,ayb,dyb)
                                         )

# BBLR-model tune shift replaced with IW model for large psix
ERR=1.e-5

def DQX(ax,ay,dx,dy,r):    
    psix=dx/ax/r
    psiy=ay/ax/r
    psiz=dy/ax
    #print("psis=",psix,psiy,psiz)
    _g   = lambda xi: g((xi/r/ax)**2,r)
    _i   = 0
    _f   = r*ax
    return -2/ax**2* integrate.quad(lambda xi:    
        xi/_g(xi)*TXL(
                 axb=xi,
                 ayb=xi*psiy/_g(xi),
                 dxb=xi*psix,
                 dyb=xi*psiz/_g(xi)), _i, _f, epsabs=ERR)[0]

def DQY(ax,ay,dx,dy,r):    
    psix=dx/ax/r
    psiy=ay/ax/r
    psiz=dy/ax
    _g   = lambda xi: g((xi/r/ax)**2,r)
    _i   = 0
    _f   = r*ax
    return  -2/ax**2/r**2* integrate.quad(lambda xi:   
    xi/_g(xi)**3*TYL(  
                axb=xi,
                ayb=xi*psiy/_g(xi),
                dxb=xi*psix,
                dyb=xi*psiz/_g(xi)), _i, _f, epsabs=ERR)[0]

# Ideal-wire real-valued FC (when one of dx, dy is zero and k is not zero)
def DmkW(m_in,k_in,ax,ay,dx,dy,r):   
   if abs(dy)<0.001:
            m=m_in
            k=k_in
            psix=abs(dx/ax/r)
            psiy=ay/ax/r
   else:
            m=k_in
            k=m_in
            psix=abs(dy/ay*r)
            psiy=ax/ay*r
   def PX(fi): return psix+sin(fi)
   def kern(fi): return (
            1/pi/k*psiy**k*
            cos((m+k)*pi/2+ m*fi)*
            (PX(fi) +sqrt(PX(fi)**2+psiy**2))**(-k)                                    )
   res=integrate.quad(kern ,0, 2*pi)[0]
   return res

# Ideal-wire real-valued FC (when one of dx, dy is zero and k=0 )
def DmkWk0(m_in,k_in,ax,ay,dx,dy,r):   
   if abs(dy)<0.001:
            m=m_in
            k=k_in
            psix=abs(dx/ax/r)
            psiy=ay/ax/r
   else:
            m=k_in
            k=m_in
            psix=abs(dy/ay*r)
            psiy=ax/ay*r
   def PX(fiy): return psix+sin(fiy)
   def kern(fiy): return (
                     1/pi/m*(sin(m*pi/2+m*fiy))*cos(fiy)*PX(fiy)
                     /( PX(fiy)**2 +psiy**2*sin(fiy)**2 ) )
   res=integrate.quad(kern ,0,2*pi)[0]
   return res

# Ideal wire tune-shifts 
def DQXW(ax,ay,dx,dy,r):
    def PX(fix): return dx+ax*r*sin(fix)
    def PY(fiy): return dy*r+ ay*sin(fiy)
    def kern(fix,fiy): return (
            -1/4/pi**2 * 2*r*sin(fix)*PX(fix)
            /( PX(fix)**2 +PY(fiy)**2 )*2/ax                           )
    return integrate.dblquad(kern ,0, 2*pi,0,2*pi)[0]

def DQYW(ax,ay,dx,dy,r):
    def PX(fix): return dx+ax*r*sin(fix)
    def PY(fiy): return dy*r+ ay*sin(fiy)
    def kern(fix,fiy): return (
            -1/4/pi**2 * 2  *sin(fiy)*PY(fiy)
            /( PX(fix)**2 +PY(fiy)**2 )*2/ay
                           )
    return integrate.dblquad(kern ,0, 2*pi,0,2*pi)[0]


# Octupole tune-shifts  
def DQXOC_FUN(ax,ay,dx,dy,r):
         return 3/2*( -2*ay**2 * r**2 +ax**2*r**4 )/dx**4
def DQYOC_FUN(ax,ay,dx,dy,r):
         return 3/2*(   ay**2  -    2*ax**2*r**2 )/dx**4

def DQXOC(ax,ay,dx,dy,r):
     if abs(dy)<0.1:
              return  DQXOC_FUN(ax,ay,dx,dy,r)
     else:
              return  DQYOC_FUN(ay,ax,dy,dx,1/r)

def DQYOC(ax,ay,dx,dy,r):
     if abs(dy)<0.1:
              return  DQYOC_FUN(ax,ay,dx,dy,r)
     else:
              return  DQXOC_FUN(ay,ax,dy,dx,1/r)



# Quadrupole tune-shifts 
def DQX0(dx,dy,r):
         return 2*r**2/dx**2

def DQY0(dx,dy,r):
         return -2/dx**2
