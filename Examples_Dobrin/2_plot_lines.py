# %%
import numpy as np
import os
dir="./"
os.listdir(dir)

# %%
dir='../Examples/Dobrin_Examples/mydata/'


import matplotlib.pyplot as plt

ipcase='ip1'
#ipcase='ip5'

pars = np.loadtxt(dir+'lrtab4py_'+ipcase+'.dat')
lrparam = np.array(pars)  
with open(dir+'names4py_'+ipcase+'.dat') as f:
        names = f.readlines()
n=len(names)
for i in range(n):
    names[i]=names[i][6:8]+names[i][10:]
    names[i]=names[i].upper()
print(names)
s=np.zeros(n)
dx=np.zeros(n)
dy=np.zeros(n)
r=np.zeros(n)
A_w_s=np.zeros(n)
B_w_s=np.zeros(n)

for i in range(n):            
        s[i],dx[i], dy[i] , r[i],A_w_s[i],B_w_s[i]=lrparam[i]


xing=145

# %%

# twofun plot
x=names

#s is for IR1 or ir5 only
# x=s

case=0#-1#0
case=-1#-1#0

dataset_1,lab1 = r,r'$r$'
dataset_2,lab2 = dx/r,r'$d_x/r=d_x^{\tt N} $'
dataset_3,lab3=dy,r'$d_y$'

if case==-1:
   ma=6
   dataset_1,lab1 = r,r'$r$'
   dataset_2,lab2 = dx/r/ma,r'$\psi_x$'
   dataset_3,lab3=dy/ma,r'$\psi_z$' 

# Creating figure 
fig = plt.figure(figsize=(14,7))

font = {'size': 12}
plt.rc('font', **font)
ax = plt.gca()
# Plotting dataset_2
#ax=fig.add_subplot(111)
plt.xticks(fontsize=11,fontweight = 'bold',rotation = 90)

plt.plot(x, dataset_2,  '-o', label=lab2,c='black', markersize=7 , markerfacecolor="none")
plt.plot(x, dataset_3, '-o', label=lab3,c='black', markersize=7, markerfacecolor="black")
#plt.grid(axis='y')

# Creating Twin axes for dataset_1
ax2 = ax.twinx()
ax2.plot(x, dataset_1, '--o', label=lab1,c='magenta', markersize=9, markerfacecolor="magenta")

# Adding title
#if case==0:
#plt.title('separatons and aspect ratios for '+str(n)+' BBLR,  '+r'$\Theta_C$ = '+str(xing)+r' $\mu$'+'rad'+"\n \n" )
#          fontweight="plain"

 
# Adding legend
ax.legend(loc='best', bbox_to_anchor=(.5, 1.2),ncol=2)
ax2.legend(loc='best', bbox_to_anchor=(.6, 1.2))
 
# adding grid
#ax.grid()
#plt.grid()

# Adding labels
if case==0:
   ax.set_xlabel(r'$\bf \leftarrow   { \tt BBLR\ names}\rightarrow$')
   ax.set_ylabel(r'$ d_x/r= \psi_x a_x$'+',    '+r'$ d_y= \psi_z  a_x$')
   ax2.set_ylabel(r'$r$')
   ax.set_ylim(-100, 100)

else:
   ax.set_xlabel(r'$\bf \leftarrow   { \tt BBLR\ names}\rightarrow$')
   ax.set_ylabel(r'$ \psi_i$')
   ax2.set_ylabel(r'$r$')
   ax.set_ylim(-5, 5)

""""
xx=n/2
down=-.5
plt.text(1/8*xx, down, 'L5_#',rotation=0, wrap=False,c='r', size=20)
plt.text(3/8*xx,down, 'R5_#',  rotation=0, wrap=False,c='g', size=20)
plt.text(5/8*xx,down, 'L1_#',  rotation=0, wrap=False,c='b', size=20)
plt.text(7/8*xx,down, 'R1_#',  rotation=0, wrap=False,c='y', size=20)
plt.text(1/4*xx,down -.0, 'IP1',rotation=0, wrap=False, size=20)
plt.text(3/4*xx,down -.0, 'IP5',  rotation=0, wrap=False, size=20)
"""
# Setting Y limits
#ax2.set_ylim(0, 35)
#ax.set_ylim(-20, 100)
plt.tight_layout()
# export plot
plt.savefig("twofun"+".pdf")


# %%


# %%
