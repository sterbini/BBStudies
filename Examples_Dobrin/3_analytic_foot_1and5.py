# %%
import numpy as np
import sys
#print(sys.path)
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

n_amp=len(coordinates['x_sig']) #total ampl points
ax_sig = coordinates['x_sig'] 
ay_sig = coordinates['y_sig'] 
plt.figure()
plt.plot(coordinates.x_sig,coordinates.y_sig,'o')
plt.axis('square')

# %%
# single call for a test
#dtune.DQXW(ax=1,ay=1,dx=1,dy=.1,r=1)


# %%
# choose the collision model  
#UseModel='BBLR'
UseModel='IW'
UseModel='OCT'
# this model defines the dtune function to be used
if UseModel=='BBLR':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQX(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQY(ax,ay,dx,dy,r)
if UseModel=='IW':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQXW(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQYW(ax,ay,dx,dy,r)
if UseModel=='OCT':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQXOC(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQYOC(ax,ay,dx,dy,r)


# %%
#USEX(ax=1,ay=1,dx=1,dy=.1,r=1)

 
# %%
import time
import itertools
colors = itertools.cycle(('r', 'g', 'b', 'y')) 
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))

MyDataDir='../Examples_Dobrin/mydata/'
for ipcase in ['ip1','ip5']:
    params = np.loadtxt(MyDataDir+'lrtab4py_'+ipcase+'.dat')
    params = np.array(params)

    with open(MyDataDir+'names4py_'+ipcase+'.dat') as f:
   	    names = f.readlines()
    n=len(names)
     
    for i in range(n):
        names[i]=names[i][6:8]+names[i][10:]
        names[i]=names[i].upper()
        
    x_tab = [[0 for col in range(n_amp)] for row in range(n)]
    y_tab = [[0 for col in range(n_amp)] for row in range(n)]
    x_tab=np.array(x_tab,dtype=np.float64)
    y_tab=np.array(y_tab,dtype=np.float64)
    
    for i in range(n):
         s,dx, dy , r,A,B=params[i]
         print("bb ",i, names[i],   " dx,dy=",dx,dy)

         s_time = time.time()
         for j in range(n_amp):
             ax = ax_sig[j] 
             ay = ay_sig[j]
             x_tab[i,j]=A**2*USEX(A*ax,B*ay,dx,dy,r)
             y_tab[i,j]=B**2*USEY(A*ax,B*ay,dx,dy,r)
             #print(ax,ay,x[i,j],y[i,j])

         e_time = time.time()
         print(f'Execution time, {(e_time-s_time):.3f} s')
    if(ipcase=='ip1'):
        x_tab1,y_tab1=x_tab,y_tab
    else:
        x_tab5,y_tab5=x_tab,y_tab
# %%
for i in range(n):
    x_individ5=x_tab5[i,:]
    y_individ5=y_tab5[i,:]
    plt.plot(
     x_individ5, y_individ5, ls=" ", lw=3, 
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
    sumx1[j]=np.sum( x_tab1[:,j])
    sumy1[j]=np.sum( y_tab1[:,j])
sumx5=np.zeros(n_amp)
sumy5=np.zeros(n_amp)
for j in range(n_amp):
    sumx5[j]=np.sum( x_tab5[:,j])
    sumy5[j]=np.sum( y_tab5[:,j])

# %%
plt.plot(
    sumx5, sumy5 , label='sum ir5',
    ls="-",marker="x"
    )
plt.title(' Model used is '+str(UseModel))
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]')
plt.grid()
plt.legend()
plt.axis('square') 

# %%

plt.plot(
    sumx1, sumy1 ,  
    sumx5, sumy5 ,  
    ls="-",marker="x")
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
plt.title('UseModel='+str(UseModel))
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]')
plt.grid()
plt.legend()
plt.axis('square') 