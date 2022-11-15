#!/usr/bin/python3
# Begin Paper Appendix 

import numpy as np
import scipy.integrate as integrate

def g(t,r):
    return np.sqrt(1+(r**2-1)*t)

def L2(u1,u2,n):
    return np.exp(-u1-u2)*integrate.quad(lambda phi:
                     np.cos(n*phi) * np.exp( u1*np.cos(phi) + u2*np.cos(2*phi)),
                                         -np.pi/2, 3/2*np.pi)[0]/2/np.pi
#def I2(u1,u2,n):
#    return integrate.quad(lambda phi:
#                     np.cos(n*phi) * np.exp( u1*np.cos(phi) + u2*np.cos(2*phi)),
#                                         -np.pi/2, 3/2*np.pi)[0]/2/np.pi
    
# kernel nominator Qmk

#def Qmz(m,azb,dzb,t):
#    u1 =  t*azb*dzb
#    u2 = -t*azb**2/4
#    return np.exp(-t*(azb**2/4+dzb**2/2))*I2(u1,u2,m)

def Qmz(m,azb,dzb,t):
    u1 =  t*azb*dzb
    u2 = -t*azb**2/4
    return np.exp(-t/2*(azb-dzb)**2)*L2(u1,u2,m)

# TS Kernel nominator in bar variables

def TSX_Bar_kern(dxb,dyb, axb,ayb,  r,t):
    return Qmz(0,ayb,dyb,t)*( Qmz(0,axb,dxb,t)+  Qmz(2,axb,dxb,t)  -2*dxb/axb* Qmz(1,axb,dxb,t)    )

def TSY_Bar_kern(dxb,dyb, axb,ayb,  r,t):
    return Qmz(0,axb,dxb,t)*( Qmz(0,ayb,dyb,t)+  Qmz(2,ayb,dyb,t)  -2*dyb/ayb* Qmz(1,ayb,dyb,t)    )

#print(TSX_Bar_kern(dxb=1,dyb=2, axb=3,ayb=4,  r=1.29203,t=1.2))

y = np.loadtxt("par.in")
dx, dy , ax , ay , r = np.array(y)
print(dx, dy , ax ,  ay , r)

########### tune shift as path integral in xi domain 

def TSX(dx,dy,ax,ay,r):    
    psix=dx/ax/r
    psiy=ay/ax/r
    print("|psix| = ",abs(psix))
    print("|psiy| = ",abs(psiy))
    return integrate.quad(lambda xi:
                          -2/ax**2*fw**2*
    TSX_Bar_kern(xi*psix, xi*dy/ax/g((xi/ax/r)**2,r),  xi,xi*psiy/g((xi/ax/r)**2,r)   ,r,1)
                          *xi/g((xi/ax/r)**2,r), 0, r*ax)[0]

def TSY(dx,dy,ax,ay,r):    
    psix=dx/ax/r
    psiy=ay/ax/r
    return integrate.quad(lambda xi:
                          -2/ax**2/r**2*fw**2*
    TSY_Bar_kern(xi*psix, xi*dy/ax/g((xi/ax/r)**2,r),  xi,xi*psiy/g((xi/ax/r)**2,r)   ,r,1)
                          *xi/g((xi/ax/r)**2,r)**3, 0, r*ax)[0]

# set fw=1 for long range, or  fw>>1 for wire
fw=1
print("fw=", fw)

if fw > 4: 
  print(" neglected flatness -- setting g(t,r)=1 when fw>>1 ")
  def g(t,r): return 1
  
print('Tune Shift X, ksi units', TSX(fw*dx,fw*dy,fw*ax,fw*ay,r)
)
print('Tune Shift Y, ksi units', TSY(fw*dx,fw*dy,fw*ax,fw*ay,r)
)


