# %%
import numpy as np
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

n_amp=len(coordinates['x_sig'])
ax_sig = coordinates['x_sig'] 
ay_sig = coordinates['y_sig'] 
plt.figure()
plt.plot(coordinates.x_sig,coordinates.y_sig,'o')
plt.axis('square')

# %%
# single call for test
dtune.DQXW(ax=1,ay=1,dx=1,dy=.1,r=1)


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
USEX(ax=1,ay=1,dx=1,dy=.1,r=1)


# %%
"""
ibb=12
nL=int(n/2)-ibb
rangen=list([nL,n-nL-1])
rangen=range(n)
"""
 
# %%
import itertools
colors = itertools.cycle(('r', 'g', 'b', 'y')) 
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))

dir='../Examples_Dobrin/mydata/'
for ipcase in ['ip1','ip5']:
#for ipcase in ['ip1']:
    pars = np.loadtxt(dir+'lrtab4py_'+ipcase+'.dat')
    lrparam = np.array(pars)

    with open(dir+'names4py_'+ipcase+'.dat') as f:
   	    names = f.readlines()
    n=len(names)
    n=1
     
    for i in range(n):
        names[i]=names[i][6:8]+names[i][10:]
        names[i]=names[i].upper()
        
    x = [[0 for col in range(n_amp)] for row in range(n)]
    y = [[0 for col in range(n_amp)] for row in range(n)]
    x=np.array(x,dtype=np.float64)
    y=np.array(y,dtype=np.float64)
    
    for i in range(n):
         s,dx, dy , r,A,B=lrparam[i]
         print("bb ",i, names[i],   " dx,dy=",dx,dy)

         s_time = time.time()
         for j in range(n_amp):
             ax = ax_sig[j] 
             ay = ay_sig[j]
             x[i,j]=A**2*USEX(A*ax,B*ay,dx,dy,r)
             y[i,j]=B**2*USEY(A*ax,B*ay,dx,dy,r)
             #print(ax,ay,x[i,j],y[i,j])

         e_time = time.time()
         print(f'Execution time, {(e_time-s_time):.3f} s')
    if(ipcase=='ip1'):
        x_ip1,y_ip1=x,y
    else:
        x_ip5,y_ip5=x,y
# %%
for i in range(n):
    xindivid5=x_ip5[i,:]
    yindivid5=y_ip5[i,:]
    plt.plot(
     xindivid5, yindivid5, ls=" ", lw=3, 
     marker='x',
     markersize=6,
      markerfacecolor=next(colors),
     markeredgecolor=next(colors),
     markeredgewidth=2,
     label=str(names[i])
      )
plt.grid()
plt.axis('square') 
    
# %%
sumx1=np.zeros(n_amp)
sumy1=np.zeros(n_amp)
for j in range(n_amp):
    sumx1[j]=np.sum( x_ip1[:,j])
    sumy1[j]=np.sum( y_ip1[:,j])
sumx5=np.zeros(n_amp)
sumy5=np.zeros(n_amp)
for j in range(n_amp):
    sumx5[j]=np.sum( x_ip5[:,j])
    sumy5[j]=np.sum( y_ip5[:,j])


# %%
        
plt.plot(
    sumx1, sumy1 , 
    ls="-",marker="x")
# %%
plt.plot(
    sumx5, sumy5 , 
    ls="-",marker="x"
    )
# %%
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]')
plt.grid()
plt.legend()
plt.axis('square') 

# %%
        
plt.plot(
    sumx1+sumx5, sumy1+sumy5, 
    ls="-",marker="x",label='ip1+ip5'
    )
plt.title('model='+str(model))
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]')
plt.grid()
plt.legend()
plt.axis('square') 


"""
    
markersize=10,
markerfacecolor="none",
markeredgecolor="red",
markeredgewidth=2,
label="analytic")
"""

    
    
    
    
    
    
    
    
    
    
        #psix=dx/r/ax
        #psiy=dy/r/ax
        #expo=1+2*psix**2+psiy**2
        #print("fact expo",expo)
        #print("around 10 psx=",dx/r/ax)
		
		
    
    
    
    
    
    
    
    
    
    
    

# %%

    
 
 
 
 
