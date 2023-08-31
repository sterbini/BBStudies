# for huge/large separations beyond coll 20 it must be zero
# %%
import matplotlib.pyplot as plt
import sys
print(sys.path)
sys.path.append('../../')
import BBStudies.Physics.Detuning_as_paper as dtune
def USEX(ax,ay,dx,dy,r):  return  dtune.DQX(ax,ay,dx,dy,r)
def USEY(ax,ay,dx,dy,r):  return  dtune.DQY(ax,ay,dx,dy,r)

# %%

#plt.style.use('seaborn-white')
import numpy as np
#from Code.code import set_global,Qmk,TX,TY,g,DQX,DQY,DQXW,DQYW,Dmk,DmkW,DQXOC,DQYOC,DQX0,DQY0

# %%
# FIG integrals
dd,tune,dd4020=True,False,False
dd,tune,dd4020=False,True,False
#dd,tune,dd4020=False,False,True

r = 1. #1.2
maxDN=80
maxa=10
#for ax=3 kernel same as g=1
gflag=0
#set_global(gflag) 
m = 4
k = 0
dy = -9
ay = .1

if tune:
    def intUSE(ax1, DN):
        return USEX(ax1,ay,DN*r,dy,r)
    word=r'$\Delta Q_x=$'
if dd:
    def intUSE(ax1,DN):
        return Dmk(m,k,ax1,ay,DN*r,dy,r)  
    word=r'$D_{40}=$'
if dd4020:
    def intUSE(ax1,DN):
        return   (Dmk(m,k,ax1,ay,DN*r,dy,r) +Dmk(2,0,ax1,ay,DN*r,dy,r))  *8/ax1**2-DQX(ax1,ay,DN*r,dy,r)
    word=''

mesh=20
axtab = np.linspace(0.001, maxa, mesh)
DNtab = -np.linspace(0.001, maxDN,3)

dataI=[]
for _DN in DNtab:
    row=[]
    for _ax in axtab:
        newval=intUSE(_ax,_DN)
        row.append(newval)
    print("=",_DN)
    print(np.max(np.abs(row)))
    plt.legend(str(_DN))
    plt.plot(row,'x')
    dataI.append(row)
plt.tight_layout()


# %%
X, Y = np.meshgrid(axtab, DNtab)
Z = dataI
from matplotlib import gridspec
gs = gridspec.GridSpec(1, 2, width_ratios=[3.5, 1]) 

CS = plt.contourf(X, Y, Z,30)


#CS = plt.contour(X, Y, Z, 40, colors='k')  # Negative contours default to dashed.

plt.clabel(CS, fontsize=14, inline=True)                                                                                    
plt.colorbar(CS)
plt.show()


# %%

CS = plt.contour(X, Y, Z, 30, colors='k')  # Negative contours default to dashed.

plt.clabel(CS, fontsize=14, inline=True)     
plt.show()


# %%
