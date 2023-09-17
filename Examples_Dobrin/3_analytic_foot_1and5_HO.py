# %%
import numpy as np
import sys
sys.path.append('../')
import BBStudies.Physics.Base as phys
# memo: here with 2D Bessel L2 (the I2 was Detuning_as_paper.py)
import BBStudies.Physics.Detuning_L2 as dtune


# %%
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sciStat
import time


# %%
coordinates = phys.polar_grid(  r_sig     = np.linspace(0.1,6.5,5),
                                theta_sig = np.linspace(0.05*np.pi/2,0.95*np.pi/2,6),
                                emitt     = [2.5e-6/7000,2.5e-6/7000])

n_amp=len(coordinates['x_sig']) #total ampl points
ax_tab = coordinates['x_sig'] 
ay_tab = coordinates['y_sig'] 
plt.figure()
plt.plot(coordinates.x_sig,coordinates.y_sig,'o')
plt.axis('square')
# %%
# Choose the collision Model (here IW is faster than BBLR)  
# it defines the dtune function USEX,Y to be used

UseModel='BBLR'
#UseModel='IW'
UseModel='OCT'

if UseModel=='BBLR':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQX(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQY(ax,ay,dx,dy,r)
if UseModel=='IW':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQXW(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQYW(ax,ay,dx,dy,r)
if UseModel=='OCT':
    def USEX(ax,ay,dx,dy,r):  return  dtune.DQXOC(ax,ay,dx,dy,r)
    def USEY(ax,ay,dx,dy,r):  return  dtune.DQYOC(ax,ay,dx,dy,r)
    
# Looping

#Flag psi_bb_as_wire serves to turn on the BBLR-as-IW approximation depending on 
# the psi parameter (inverted distance to the strong beam core)
# Set this flag to a huge value to ignore. 

#psi_bb_as_wire=1.e100



MyDataDir='../Examples_Dobrin/mydata/'
#MyDataDir='../Examples_Dobrin/mydata_NOS/'
i_bb=0

#for ipcase in ['ip1','ip5']:
#for ipcase in ['ip1ho','ip5ho']:
for ipcase in [ 'ip1','ip5','ip1ho','ip5ho']:
    params = np.loadtxt(MyDataDir+'data_'+ipcase+'.dat')
 
    params = np.array(params)
    with open(MyDataDir+'names_'+ipcase+'.dat') as f:
   	    names = f.readlines()
    nbb=len(names)
    split_factor=1
    if 'ho' in ipcase:
        split_factor=1/nbb
        
        def USEX(ax,ay,dx,dy,r):  return  dtune.DQX(ax,ay,dx,dy,r)
        def USEY(ax,ay,dx,dy,r):  return  dtune.DQY(ax,ay,dx,dy,r)
    for i in range(nbb):
        names[i]=names[i][6:8]+names[i][10:13]
        names[i]=names[i].upper()
        
    x_tab = [[0 for jj in range(n_amp)] for ii in range(nbb)]
    y_tab = [[0 for jj in range(n_amp)] for ii in range(nbb)]
    x_tab=np.array(x_tab,dtype=np.float64)
    y_tab=np.array(y_tab,dtype=np.float64)
    
    for i in range(nbb):
         i_bb=i_bb+1
         s,dx, dy , r,A,B=params[i]
         print(ipcase,"i_bb ",i_bb, names[i],   " dx,dy=",dx,dy)
         

         s_time = time.time()
         for j in range(n_amp):
             ax = ax_tab[j] 
             ay = ay_tab[j]
             #psix_test=np.abs(dx/ax/r)
             #psiz_test=np.abs(dy/ax)
             
             resx=A**2*USEX(A*ax,B*ay,dx,dy,r)
             resy=B**2*USEY(A*ax,B*ay,dx,dy,r)
             #if psix_test > psi_bb_as_wire:
             #    resx=A**2*dtune.DQXW(A*ax,B*ay,dx,dy,r)
             #    resy=B**2*dtune.DQYW(A*ax,B*ay,dx,dy,r)            
             #else:
             #    resx=A**2*USEX(A*ax,B*ay,dx,dy,r)
             #    resy=B**2*USEY(A*ax,B*ay,dx,dy,r)


             x_tab[i,j]=split_factor*resx
             y_tab[i,j]=split_factor*resy

             #print(ax,ay,x_tab[i,j],y_tab[i,j])

         e_time = time.time()
         print(f'Execution time, {(e_time-s_time):.3f} s')
        
    if(ipcase=='ip1'):
        names1,x_tab1,y_tab1=names,x_tab,y_tab
    if(ipcase=='ip5'):
        names5,x_tab5,y_tab5=names,x_tab,y_tab
    if(ipcase=='ip1ho'):
        names1ho,x_tab1ho,y_tab1ho=names,x_tab,y_tab
    if(ipcase=='ip5ho'):
        names5ho,x_tab5ho,y_tab5ho=names,x_tab,y_tab
 # %%
 
 #
 # function to plot each 
 #
def plot_indiv(names,ccc,xtab,ytab):
    for i in range(len(names)):
        if "L" in names[i]: 
                cc=ccc[0]
        if "R" in names[i]: 
                cc=ccc[1]
        x_individ=xtab[i,:]
        y_individ=ytab[i,:]
        plt.plot(x_individ, y_individ, ls="--", lw=1, 
            marker='x',markersize=3,color=cc  
        )
    #plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))
    plt.grid()
    plt.axis('square')
    sumx=np.zeros(n_amp)
    sumy=np.zeros(n_amp)
    for j in range(n_amp):
        sumx[j]=np.sum( xtab[:,j])
        sumy[j]=np.sum( ytab[:,j])
    #plt.plot(
    #sumx, sumy , label='sum',
    #ls="-",marker="x",color=cc  )
    return sumx,sumy
    
    
# %%
##########
#  PLOTS
##########
colors = (('r', 'g'), ('b','yellow'))

sumx5, sumy5=plot_indiv(names5,colors[0],x_tab5,y_tab5)
sumx1, sumy1=plot_indiv(names1,colors[1],x_tab1,y_tab1)
sumx1ho, sumy1ho=plot_indiv(names1ho,colors[1],x_tab1ho,y_tab1ho)
sumx5ho, sumy5ho=plot_indiv(names5ho,colors[0],x_tab5ho,y_tab5ho)
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]') 
plt.grid()
plt.legend()
plt.axis('square')

# %%
plt.plot(sumx1+sumx5+sumx1ho+sumx5ho, sumy1+sumy5+sumy1ho+sumy5ho,ls="--",marker="x",
         label='ip1+ip5', markersize=4,color='black' )
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]') 
plt.grid()
plt.legend()
plt.axis('square')
plt.title(' Model used is '+str(UseModel))





# %%
#sumx1, sumy1=plot_indiv(names1,colors[1],x_tab1,y_tab1)

# %%
#sumx5, sumy5=plot_indiv(names5,colors[0],x_tab5,y_tab5)

# %%
#sumx1ho, sumy1ho=plot_indiv(names1ho,colors[1],x_tab1ho,y_tab1ho)
# %%
#sumx5ho, sumy5ho=plot_indiv(names5ho,colors[0],x_tab5ho,y_tab5ho)



# %%
