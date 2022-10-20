
"""
Created on Tue May 17 2022
@author: pbelange
Description :
"""

import numpy as np
import scipy.integrate as integrate
import scipy.special as sciSpec




#================================================================================
#    Tune Shift
#================================================================================
def DQx_DQy(ax,ay,dx_sig,dy_sig,A_w_s,B_w_s,r,xi,fw=1):
    """
    Notes: 
    The function expects an array for ax,ay, and a single value for the other parameters
    --------
    ax,ay         -> normalized amplitude, ax = x/sigma_weak
    r             -> sigma_y/sigma_x
    dx_sig,dy_sig -> normalized bb separation, dx_sig = dx/sigma_strong
    xi            -> beam-beam parameter
    A_w_s         -> sigma_w_x/sigma_s_y
    B_w_s         -> sigma_w_y/sigma_s_x
    fw            -> reduction factor, sig_x,y -> sig_x,y/fw
    """
    ax = np.array(ax)
    ay = np.array(ay)

    DQx = xi*((A_w_s*fw)**2)*np.array([DQx_norm(_ax*(A_w_s*fw),_ay*(B_w_s*fw),dx_sig*fw,dy_sig*fw,r) for _ax,_ay in zip(ax,ay)])
    DQy = xi*((B_w_s*fw)**2)*np.array([DQy_norm(_ax*(A_w_s*fw),_ay*(B_w_s*fw),dx_sig*fw,dy_sig*fw,r) for _ax,_ay in zip(ax,ay)])
    return DQx,DQy
#================================================================================



#================================================================================
#    Bessel 2D 
#================================================================================

def Bess2D_kern(phi,u1,u2,n):
    return np.cos(n*phi) * np.exp( u1*np.cos(phi) + u2*np.cos(2*phi))

def Bess2D(u1,u2,n):
    _i = -np.pi/2
    _f = 3/2*np.pi
    return np.exp(-u1-u2)/2/np.pi*integrate.quad(lambda phi: 
                                    Bess2D_kern(phi,u1,u2,n),_i,_f)[0]
#================================================================================


#================================================================================
#    Kernel for the path integral 
#================================================================================
def Qmz(azb,dzb,n,t=1):
    u1 =  t*azb*dzb
    u2 = -t*azb**2/4
    return np.exp(-t/2*(azb-dzb)**2)*Bess2D(u1,u2,n)
#================================================================================



#================================================================================
#    Path integral for tune shift
#================================================================================

def g(t,r):
    return np.sqrt(1+(r**2-1)*t)

def TSX_kern(axb,ayb,dxb,dyb):
    return Qmz(ayb,dyb,0)*( Qmz(axb,dxb,0) + Qmz(axb,dxb,2) - 2*dxb/axb* Qmz(axb,dxb,1) )

def TSY_kern(axb,ayb,dxb,dyb):
    return Qmz(axb,dxb,0)*( Qmz(ayb,dyb,0) + Qmz(ayb,dyb,2) - 2*dyb/ayb* Qmz(ayb,dyb,1) )


#--------------------------------------------------------------------------------
def DQx_norm(ax,ay,dx,dy,r): 
    """The function expects a single value for all parameters"""   
    psix = dx/ax/r
    psiy = ay/ax/r
    psiz = dy/ax
    
    _g   = lambda phi: g((phi/r/ax)**2,r)

    _i   = 0
    _f   = r*ax
    return (-2/ax**2) * integrate.quad(lambda phi:
                                    phi*(1/_g(phi))*TSX_kern(   axb = phi,
                                                                ayb = phi*psiy/_g(phi),
                                                                dxb = phi*psix,
                                                                dyb = phi*psiz/_g(phi)),_i,_f)[0]
#--------------------------------------------------------------------------------
def DQy_norm(ax,ay,dx,dy,r):
    """The function expects a single value for all parameters"""     
    psix = dx/ax/r
    psiy = ay/ax/r
    psiz = dy/ax
    
    _g   = lambda phi: g((phi/r/ax)**2,r)

    _i   = 0
    _f   = r*ax
    return (-2/ax**2) * integrate.quad(lambda phi:
                                    phi*(1/r**2/_g(phi)**3)*TSY_kern(   axb = phi,
                                                                        ayb = phi*psiy/_g(phi),
                                                                        dxb = phi*psix,
                                                                        dyb = phi*psiz/_g(phi)),_i,_f)[0]
#================================================================================





#================================================================================
#    Octupolar approximation
#================================================================================

def BBLR_octupole(Jx,Jy,betx,bety,k1,k3):
    
    # Quadrupole contribution
    DQx =  1/(4*np.pi)*k1*betx
    DQy = -1/(4*np.pi)*k1*bety
        
    # Octupole contribution
    DQx += 3/(8*np.pi)*(k3/np.math.factorial(3))*(betx**2 * Jx - 2*betx*bety*Jy)
    DQy += 3/(8*np.pi)*(k3/np.math.factorial(3))*(bety**2 * Jy - 2*bety*betx*Jx)
    
    return DQx,DQy








def Z1(x):
    return np.exp(-x)*(sciSpec.iv(0,x)-sciSpec.iv(1,x))

def Z2(x):
    return np.exp(-x)*sciSpec.iv(0,x)

def HeadOn_round_generating(t,alphax,alphay,r):
    # Assume in x direction. For y, change x->y, y->x, r->1/r
    # See: https://inspirehep.net/files/d7fed02b4b59558edf043a65f0f92049
    prefactor = 1/((1+t)**(3/2)* (1+t/r**2)**(1/2))
    return (1+1/r)/2* prefactor * Z1(alphax/(1+t)) * Z2(alphay/(1+t/r**2))
    
    
def HeadOn_round(ax,ay,r,emitt,xi):
    
    if isinstance(emitt, (int, float)):
        emitt = emitt*np.ones(2)
    if isinstance(xi, (int, float)):
        xi = xi*np.ones(2)
        
    # Normalization from Chao
    alphax = ax**2/(2*(1+emitt[1]/emitt[0]))
    alphay = ay**2/(2*(1+emitt[0]/emitt[1]))
    
    # Tune shifts, t is the integration variable
    DQx_n = np.array([integrate.quad(lambda t: HeadOn_round_generating(t,_alphax,_alphay,r)  , 0, np.inf)[0] for _alphax,_alphay in zip(alphax,alphay)])
    DQy_n = np.array([integrate.quad(lambda t: HeadOn_round_generating(t,_alphay,_alphax,1/r), 0, np.inf)[0] for _alphax,_alphay in zip(alphax,alphay)])
    
    return -xi[0]*DQx_n,-xi[1]*DQy_n


#================================================================================
#    LEGACY DETUNING
#================================================================================
class legacy():

    def __init():
        pass
    @classmethod
    def DQx_DQy(self,ax,ay,r,dx_n,dy_n,xi,A_w_s,B_w_s,method='fast'):
        """
        Notes: 
        The function expects an array for ax,ay, and a single value for the other parameters
        --------
        ax,ay -> normalized amplitude, ax = x/sigma_weak
        r     -> sigma_y/sigma_x
        dx,sy -> normalized bb separation, dx_n = dx/sigma_strong
        xi    -> beam-beam parameter
        A_w_s -> sigma_w_x/sigma_s_y
        B_w_s -> sigma_w_y/sigma_s_x
        """
        DQx = np.array([2*xi*self.dC00dx(A_w_s*_ax,B_w_s*_ay,r,dx_n,dy_n,method)/(A_w_s*_ax) for _ax,_ay in zip(ax,ay)])
        DQy = np.array([2*xi*self.dC00dy(A_w_s*_ax,B_w_s*_ay,r,dx_n,dy_n,method)/(B_w_s*_ay) for _ax,_ay in zip(ax,ay)])
        return (A_w_s**2)*DQx,(B_w_s**2)*DQy

    #================================================================================
    #---------------------------------------
    @classmethod
    def g(self,t,r):
        return np.sqrt(1+(r**2-1)*t)
    #---------------------------------------
    #---------------------------------------
    @classmethod
    def Bess2D_generating_real(self,phi,X,Y,n):
        arg =  - X*np.sin(phi) + 2*Y*np.sin(phi)**2
        return np.cos(n*phi)*np.exp(arg)
    @classmethod
    def Bess2D_generating_imag(self,phi,X,Y,n):
        arg =  - X*np.sin(phi) + 2*Y*np.sin(phi)**2
        return -np.sin(n*phi)*np.exp(arg)

    @classmethod
    def Bess2D_INT(self,X,Y,n):

        # Choosing between the real part or the imaginary part depending on n
        generatingFun,sign = {0:(self.Bess2D_generating_real,(-1)**(n/2)),
                            1:(self.Bess2D_generating_imag,(-1)**((n-1)/2))}[n%2]
        # computing the integral
        integratedFactor = integrate.quad(lambda phi: generatingFun(phi,X,Y,n), 0, 2*np.pi)[0]

        return sign*np.exp(-X-2*Y)/2/np.pi * integratedFactor

    @classmethod
    def Bess2D_SUM(self,X,Y,n):
        qmax = 40
        order = np.arange(-qmax,qmax+1)
        
        q,XX = np.meshgrid(order,X)
        _,YY = np.meshgrid(order,Y)
        
        return np.exp(-X-Y)*np.sum(sciSpec.iv(n-2*q,XX)*sciSpec.iv(q,YY),axis=1)

    @classmethod
    def Bess2D(self,X,Y,n,method='int'):
        return {'int':Bess2D_INT,'sum':Bess2D_SUM}[method](X,Y,n)





    #---------------------------------------
    #---------------------------------------
    @classmethod
    def Q0z(self,t,azbar,dzbar):
        X =  t*azbar*dzbar
        Y = -t*azbar**2/4
        return np.exp(-t/2*(azbar-dzbar)**2)*self.Bess2D(X,Y,0)
    #---------------------------------------
    #---------------------------------------
    @classmethod
    def dQ0daz(self,t,azbar,dzbar,etaz):
        X    =  t*azbar*dzbar
        Y    = -t*azbar**2/4
        return etaz*np.exp(-t/2*(azbar-dzbar)**2)*(-azbar/2*(self.Bess2D(X,Y,0)+self.Bess2D(X,Y,2)) + dzbar*self.Bess2D(X,Y,1))
    #---------------------------------------
    #---------------------------------------
    @classmethod
    def dQ0daz_Bess_generating(self,phi,X,Y,azbar,dzbar):
        arg =  -X*np.sin(phi) + 2*Y*np.sin(phi)**2
        # exp(arg) is commong to Bess_0, Bess_1, Bess_2. The prefactor coming from the sum is:
        pre = azbar*(np.cos(2*phi)-1)/2 - dzbar*np.sin(phi)
        return pre*np.exp(arg)
    @classmethod
    def Fast_dQ0daz(self,t,azbar,dzbar,etaz):
        X    =  t*azbar*dzbar
        Y    = -t*azbar**2/4
        
        integratedFactor = integrate.quad(lambda phi: self.dQ0daz_Bess_generating(phi,X,Y,azbar,dzbar), 0, 2*np.pi)[0]
        
        return etaz*np.exp(-t/2*(azbar-dzbar)**2)*(np.exp(-X-2*Y)/2/np.pi)*integratedFactor
    #---------------------------------------


    #---------------------------------------
    @classmethod
    def dC00dx_generating(self,t,ax,ay,r,dx_n,dy_n,method='regular'):
        # bar variables
        axbar,aybar,dxbar,dybar = ax*r, ay/g(t,r) , dx_n , r*dy_n/g(t,r)

        # modified version of eta to include all the prefactors
        etax_modif = r/g(t,r)
        
        derivative = {'regular':self.dQ0daz,'fast':self.Fast_dQ0daz}[method]

        return derivative(t,axbar,dxbar,etax_modif)*self.Q0z(t,aybar,dybar)

    @classmethod
    def dC00dy_generating(self,t,ax,ay,r,dx_n,dy_n,method='regular'):
        # bar variables
        axbar,aybar,dxbar,dybar = ax*r, ay/g(t,r) , dx_n , r*dy_n/g(t,r)

        # modified version of eta to include all the prefactors
        etay_modif = 1/(g(t,r)**2)
        
        derivative = {'regular':self.dQ0daz,'fast':self.Fast_dQ0daz}[method]

        return derivative(t,aybar,dybar,etay_modif)*Q0z(t,axbar,dxbar)

    @classmethod
    def dC00dx(self,ax,ay,r,dx_n,dy_n,method = 'regular'):
        return integrate.quad(lambda t: self.dC00dx_generating(t,ax,ay,r,dx_n,dy_n,method), 0, 1)[0]

    @classmethod
    def dC00dy(self,ax,ay,r,dx_n,dy_n,method = 'regular'):
        return integrate.quad(lambda t: self.dC00dy_generating(t,ax,ay,r,dx_n,dy_n,method), 0, 1)[0]
    #---------------------------------------
    #================================================================================
    #================================================================================



