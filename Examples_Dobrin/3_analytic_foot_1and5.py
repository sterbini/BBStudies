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
coordinates = phys.polar_grid(  r_sig     = np.linspace(1,6.5,5),
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
UseModel='IW'
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
    
# %%
psi_bb_as_wire=1.e100

MyDataDir='../Tests_Dobrin/mydata/'
i_bb=0

for ipcase in ['ip1','ip5']:
    params = np.loadtxt(MyDataDir+'data_'+ipcase+'.dat')
#    params = np.loadtxt(MyDataDir+'data_'+ipcase+'_NOS.dat')

    params = np.array(params)
    with open(MyDataDir+'names_'+ipcase+'.dat') as f:
   	    names = f.readlines()
    nbb=len(names)
    n_bb=nbb
    for i in range(nbb):
        names[i]=names[i][6:8]+names[i][10:13]
        names[i]=names[i].upper()
        
    x_tab = [[0 for col in range(n_amp)] for row in range(nbb)]
    y_tab = [[0 for col in range(n_amp)] for row in range(nbb)]
    x_tab=np.array(x_tab,dtype=np.float64)
    y_tab=np.array(y_tab,dtype=np.float64)
    
    for i in range(nbb):
         i_bb=i_bb+1
         s,dx, dy , r,A,B=params[i]
         print("i_bb ",i_bb, names[i],   " dx,dy=",dx,dy)
         

         s_time = time.time()
         for j in range(n_amp):
             ax = ax_tab[j] 
             ay = ay_tab[j]
             psix_test=np.abs(dx/ax/r)
             #psiz_test=np.abs(dy/ax)
             
             if psix_test > psi_bb_as_wire:
                 resx=A**2*dtune.DQXW(A*ax,B*ay,dx,dy,r)
                 resy=B**2*dtune.DQYW(A*ax,B*ay,dx,dy,r)            
             else:
                 resx=A**2*USEX(A*ax,B*ay,dx,dy,r)
                 resy=B**2*USEY(A*ax,B*ay,dx,dy,r)

#             resx=A**2*USEX(A*ax,B*ay,dx,dy,r)
#             resy=B**2*USEY(A*ax,B*ay,dx,dy,r)
#             resxIW=A**2*dtune.DQXW(A*ax,B*ay,dx,dy,r)
#             resyIW=B**2*dtune.DQYW(A*ax,B*ay,dx,dy,r)
#             errx_BBasIW=np.abs((resxIW-resx)/resx)
#             erry_BBasIW=np.abs((resyIW-resy)/resy)
#             #must be 2 for 10
#             if errx_BBasIW > 0.001  and psix_test > 10 and ax >2:
#                 print(f'large psi but NOT as IW  |psix|=, {psix_test:.3f} ')
#                 print(f'                                          |psiz|=, {psiz_test:.3f} ')
#                 print(ax,ay)
#                 print(f'errx =, {  errx_BBasIW :.3e} ')
#                 print(f'erry =, {  erry_BBasIW :.3e} ')

             x_tab[i,j]=resx
             y_tab[i,j]=resy

#             print(ax,ay,x_tab[i,j],y_tab[i,j])

         e_time = time.time()
         print(f'Execution time, {(e_time-s_time):.3f} s')
    if(ipcase=='ip1'):
        names1,x_tab1,y_tab1=names,x_tab,y_tab
    else:
        names5,x_tab5,y_tab5=names,x_tab,y_tab
# %%

colors = ('r', 'g', 'b','y')
import itertools
for i in range(nbb):
    if "L5" in names5[i]: 
            cc=colors[0]
    if "R5" in names5[i]: 
            cc=colors[1]
    if "L1" in names5[i]: 
            cc=colors[2]
    if "R1" in names5[i]: 
            cc=colors[3]
    x_individ5=x_tab5[i,:]
    y_individ5=y_tab5[i,:]
    plt.plot(x_individ5, y_individ5, ls="-", lw=1, 
        marker='x',markersize=3,color=cc  
      )
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))
plt.grid()
plt.axis('square')     
plt.show()

# %%
for i in range(nbb):
    if "L5" in names1[i]: 
            cc=colors[0]
    if "R5" in names1[i]: 
            cc=colors[1]
    if "L1" in names1[i]: 
            cc=colors[2]
    if "R1" in names1[i]: 
            cc=colors[3]   
    x_individ1=x_tab1[i,:]
    y_individ1=y_tab1[i,:]
    plt.plot(x_individ1, y_individ1, ls="-", lw=1, 
     marker='o',markersize=3 ,color=cc
      
      )
plt.ticklabel_format(style='sci', axis='both', scilimits=(-2,-2))
plt.grid()
plt.axis('square') 
plt.show()
print(names1)
print(names5)

    
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
# %%

plt.show()

# %%
