# %%
import numpy as np
dir='/home/kaltchev/BBStudies_old/Examples/AsPaper/mydata/'
#import os
#os.listdir(dir)

# %%
import sys
print(sys.path)
sys.path.append('../../')
import BBStudies.Physics.Base as phys
import BBStudies.Physics.Detuning_as_paper as dtune

# %%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sciStat


# %%

coordinates = phys.polar_grid(  r_sig     = np.linspace(1,6.5,5),
                                theta_sig = np.linspace(0.05*np.pi/2,0.95*np.pi/2,6),
                                emitt     = [2.5e-6/7000,2.5e-6/7000])

npts=len(coordinates['x_sig'])
ax_sig = coordinates['x_sig'] 
ay_sig = coordinates['y_sig'] 
plt.figure()
plt.plot(coordinates.x_sig,coordinates.y_sig,'o')
plt.axis('square')
# %%

# %%
dtune.DQX(1,1,1,1,1)
# %%
import time
model='BBLR'
model='IW'
model1='OCT'

if model=='BBLR':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQX(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQY(ax,ay,dx,dy,r)
if model=='IW':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQXW(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQYW(ax,ay,dx,dy,r)
if model=='OCT':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQXOC(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQYOC(ax,ay,dx,dy,r)

# %%
ipcase='ip1'
ipcase='ip5'
pars = np.loadtxt(dir+'lrtab4py_'+ipcase+'.dat')
lrparam = np.array(pars)  
with open(dir+'names4py_'+ipcase+'.dat') as f:
    names = f.readlines()
n=len(names)
for i in range(n):
    names[i]=names[i][6:8]+names[i][10:]
    names[i]=names[i].upper()
print(names)


# %%
x = [[0 for col in range(npts)] for row in range(n)]
y = [[0 for col in range(npts)] for row in range(n)]
x=np.array(x,dtype=np.float64)
y=np.array(y,dtype=np.float64)

############s
for i in range(n):    
#for i in range(1):    

    s,dx, dy , r,A,B=lrparam[i]
    print("bb ",i, names[i], " of ", n, "dx dy=",dx,dy)
    s_time = time.time()
    for j in range(npts):
        ax = ax_sig[j] #double to see
        ay = ay_sig[j]
        x[i,j]=A**2*USEX(A*ax,B*ay,dx,dy,r)
        y[i,j]=B**2*USEY(A*ax,B*ay,dx,dy,r)
        #print(ax,ay,x[i,j],y[i,j])
        #psix=dx/r/ax
        #psiy=dy/r/ax
        #expo=1+2*psix**2+psiy**2
        #print("fact expo",expo)
        #print("around 10 psx=",dx/r/ax)
    e_time = time.time()
    print(f'Execution time, {(e_time-s_time):.3f} s')


# %%
sumx=np.zeros(npts)
sumy=np.zeros(npts)
for j in range(npts):
  sumx[j]=np.sum( x[:,j])
  sumy[j]=np.sum( y[:,j])
   

# %%

import itertools
colors = itertools.cycle(('r', 'g', 'b', 'y')) 
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))

for i in range(n):    
    xbbn=x[i,0:npts]
    ybbn=y[i,0:npts]
    plt.plot(
        xbbn, ybbn, ls=" ", lw=3, 
        marker="x",
        markersize=6,
        markerfacecolor=next(colors),
#        markeredgecolor=next(colors),
        markeredgewidth=2,
        label=str(names[i])
          )
 

# %%
 
plt.title('model='+str(model)+ ' reset='+str(0))
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]')
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))
 
plt.plot(
    sumx, sumy , ls="-",
    marker="x",
    markersize=10,
    markerfacecolor="none",
    markeredgecolor="red",
    markeredgewidth=2,
    label="analytic"
)
plt.grid()
plt.axis('square')


# %%
