# %%
import numpy as np
dir='/home/kaltchev/BBStudies_old/Examples/AsPaper/mydata/'
# %%
ipcase='ip1'
#ipcase='ip5'
# %%
import sys
print(sys.path)
sys.path.append('../../')
import BBStudies.Physics.Base as phys
import BBStudies.Physics.Detuning_as_paper as dtune
# %%
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
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sciStat


# %%
""""
coordinates = phys.polar_grid(  r_sig     = np.linspace(1,6.5,5),
                                theta_sig = np.linspace(0.05*np.pi/2,0.95*np.pi/2,6),
                                emitt     = [2.5e-6/7000,2.5e-6/7000])

namp=len(coordinates['x_sig'])
ax_sig = coordinates['x_sig'] 
ay_sig = coordinates['y_sig'] 
plt.figure()
plt.plot(coordinates.x_sig,coordinates.y_sig,'o')
plt.axis('square')
"""
# %%
ax_sig=[6,.1]
ay_sig=[.1,6]
ax_sig=[6]
ay_sig=[.1]

namp=len(ax_sig)
# %%
qx = [[0 for col in range(namp)] for row in range(n)]
qy = [[0 for col in range(namp)] for row in range(n)]
qx=np.array(qx,dtype=np.float64)
qy=np.array(qy,dtype=np.float64)

# %%
dtune.DQX(1,1,1,1,1)
# %%
import time
model='BBLR'
model1='IW'
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

############
for i in range(n):    

    s,dx, dy , r,A,B=lrparam[i]
    A=1
    B=1
    print("bb ",i, names[i], " of ", n, "dx dy=",dx,dy)
    s_time = time.time()
    for j in range(namp):
        ax = ax_sig[j] #double to see
        ay = ay_sig[j]
        qx[i,j]=A**2*(USEX(A*ax,B*ay,dx,dy,r))-0*USEX(.0001,.0001,dx,dy,r)
        qy[i,j]=A**2*(USEY(A*ax,B*ay,dx,dy,r))-0*USEY(.0001,.0001,dx,dy,r)
    e_time = time.time()
    print(f'Execution time, {(e_time-s_time):.3f} s')


# %%
import itertools
colors = itertools.cycle(('g', 'r', 'b','y','black')) 
fig = plt.figure(figsize=(7,5))
#plt.xticks(fontsize=11,rotation = 90)
#plt.xticks(np.arange(1, 19, 17))
#plt.yticks(fontsize=11)
axL=fig.add_subplot(111)

axR = axL.twinx()
Ltxt='BBLR MEF'+r'$_x$'
xing1="145"
cc=next(colors)
axL.legend(loc='lower left',bbox_to_anchor=(0, 1.1),ncol=1,
           title='Contrib. along '+ipcase,  borderaxespad=0, frameon=False)
axR.legend(loc='lower left',bbox_to_anchor=(.7, 1.1),ncol=1,
           title=r'$a_x=$'+ str(ax),    borderaxespad=0, frameon=False)
x= list(range(1,n+1))
axL.set_xlabel(r'$ j $',size=20)
axR.set_xlabel(r'$ j $',size=20)
for j in range(namp):
    cc=next(colors)
    axL.plot(x, qx[:,j],'o',  ls="-",lw=3   , 
             label='Qx '+'ax='+str(ax_sig[j])+str(ay_sig[j]),c=cc, markersize=11)
    axR.plot(x, qy[:,j], '+', ls="--",lw=3 ,
               label='Qy '+'ax='+str(ax_sig[j]),c=cc, markersize=11)
axL.set_xticks([1,7,n])
      
plt.legend()
plt.grid()
plt.show()


     
# %%