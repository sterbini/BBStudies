# %%

#########################
#  
# mini-ring footprint in two ways:
# analyitic and with get_footprint()
# 
#  The "mini" ring lattice is
#  made in MadX w/o beam-beam, and saved as a thin sequence. 

#  Thus the weak beam, the xing bump and the IP locaton are fixed. 

# A sinlge BBLR is installed at a distance s_bb from the IP.

# If flag HO_case=False, then the Strong-beam orbit x_str and the collision 
# is chosen arbitrary. The full separation is sepx=x_wk-x_str.

# If flag HO_case=True, then s_bb = 0, x_str=0, so it becomes Head-On

# The strong beam rms at the collision is next defined for 3 cases
#    case=0  exactly assymetric
#    case=1 arbitrary rms
#    case =2 to get a desired B_w_s and flatness r 

###########################

import xtrack as xt
from cpymad.madx import Madx
import sys
sys.path.append('../')

#import xtrack as xt
import xpart as xp
import xobjects as xo
import numpy as np
from numpy import *
import matplotlib.pyplot as plt

import BBStudies.Physics.Base as phys
#import BBStudies.Plotting.BBPlots as bbplt
#import BBStudies.Physics.Constants as cst
#import xpart as xp
import xfields as xf

import BBStudies.Physics.Detuning_L2 as dtuneL2
from mad.plot_twiss import *

#  
energy = 1000e9  
#

# %%
#
# and do the mad run
#
mad = Madx()
mad.call("mad/myscript.madx")
mad.call("mad/call_mini_seq.madx")
mad.call("mad/seqb1thin")
mad.use("lhcb1")


# Build Xtrack line importing MAD-X expressions
mini = xt.Line.from_madx_sequence(mad.sequence['lhcb1'], deferred_expressions=True)

s_ip=mini.get_s_position('ip')

######################### 
#
#  Begin lattice and bb setup
#  
########################### 
#
#   set the bump 
#   should match it  

mini.vars['kxl1b1']= -1.41220e-03   
mini.vars['kxl2b1']= -6.06038e-05
mini.vars['kxr1b1']=  6.61748e-04
mini.vars['kxr2b1']= -8.29579e-04


# %%
HO_case=False   # If true, will install BBLR at the IP so is is HO

s_bb=-4  # BBLR distance from the IP

if HO_case:
    s_bb=0

#mini.particle_ref = xp.Particles(
#                    mass0=xp.PROTON_MASS_EV, q0=1, energy0=energy)

context = xo.ContextCpu()         # For CPU
#
# Define a dummy bb element - 
# there may be some better way
#
bb_any = xf.BeamBeamBiGaussian2D(
        _context=context,
#	    other_beam_num_particles=123,
        other_beam_q0 = 1,
        other_beam_shift_x=.123,            
        ref_shift_x=0,
        other_beam_beta0 = [123,123],
        other_beam_Sigma_11 = .123,
        other_beam_Sigma_33 = .123,
        other_beam_num_particles=0)

# Install the BBLR 

bbname='bbkick'

mini.insert_element(element=bb_any,name=bbname, at_s=s_ip+s_bb)
names=mini.element_names

mini.build_tracker()

# %%
rp = 1.534698E-18
########
Nb=1.e13
x_str=0.004     #  must be > 0

bunch_intensity = Nb 
nemitt=6e-5
nemitt_x = nemitt 
nemitt_y = nemitt# 
gamma = energy/xp.PROTON_MASS_EV
physemit_x = nemitt_x/gamma
physemit_y = nemitt_y/gamma
ksi= (bunch_intensity*rp)/(4*np.pi*nemitt_x)
print('ksi=',{ksi})
# %%

mini.particle_ref = xp.Particles(
    mass0=xp.PROTON_MASS_EV, q0=1, energy0=energy )

tw = mini.twiss(method='4d')

# Initial tune and chromaticity values
print(f"Qx = {tw['qx']:.5f} Qy = {tw['qy']:.5f} "
      f"Q'x = {tw['dqx']:.5f} Q'y = {tw['dqy']:.5f}")

pltwis(tw)
# %%

#### Take strong from weak ###
names=mini.element_names

for ind,ee in enumerate(names):
	if ee.find(bbname)!=-1:
		betx_wk=tw['betx'][ind]
		bety_wk=tw['bety'][ind]
		x_wk=tw['x'][ind] 
		print(mini[ee])

if HO_case:
    x_str=0

sepx=x_wk-x_str
sigx_wk=np.sqrt(physemit_x *betx_wk)#
sigy_wk=np.sqrt(physemit_y *bety_wk)#

# choose strong beam rms either exactly 
# asymmetrcic or arbitrary

case=2
if HO_case:
    case=0
    
if case==0:
    sigx_str=sigy_wk 
    sigy_str=sigx_wk 

if case==1:
# Choose arbitrary rms for weak (weak rms are already fixed)
    sigx_str=0.00123
    sigy_str=0.00234

if case==2:
# Alternatively choose seom desired r and B_w_s
    B_want=1.234
    r_want=1.23
    sigx_str=1/B_want*sigy_wk
    sigy_str=r_want*sigx_str


for my_index,ee in enumerate(names):
	if ee.find(bbname)!=-1:
            mini[ee].other_beam_q0 = 1,
            mini[ee].other_beam_shift_x=x_str,
            mini[ee].ref_shift_x=0,
            mini[ee].other_beam_beta0 = [123,123]  
            mini[ee].other_beam_Sigma_11 = sigx_str**2 , 
            mini[ee].other_beam_Sigma_33 = sigy_str**2 ,
            mini[ee].other_beam_num_particles=Nb


for my_index,ee in enumerate(names):
	if ee.find(bbname)!=-1:
              print(mini[ee])



Qx_0=tw['qx']    #tw['mux'][-1]
Qy_0=tw['qy']    #tw['muy'][-1]

# check that this is small
#twbb = mini.twiss(method='4d')
#print(twbb['qx']-tw['qx'])
        
#Qx_bb=twbb['qx']    #twbb['mux'][-1]
#Qy_bb=twbb['qy']    #twbb['muy'][-1]

######################### 
#
#  End of lattice and bb setup
#  
########################### 

 
# %%

# %%
# Generate ax-ay set 
# hope it is the same as in get_footprint
a_min = .05
a_max = 6
n_amp = 10
n_wings = 10
coordinates = phys.polar_grid(r_sig=np.linspace(a_min, a_max, n_amp),
                              theta_sig=np.linspace(
                                  0.05*np.pi/2, 0.95*np.pi/2, n_wings),
                              emitt=[nemitt/energy, nemitt/energy])
ax = coordinates['x_sig']
ay = coordinates['y_sig']


# %%
#
# Analytic detuning with Model=BBLR
#
from multiprocessing import Pool
Model='BBLR'

fw = 1
dx_sig=sepx/sigx_str
dy_sig=0
r = sigy_str/sigx_str
A_w_s = sigx_wk/sigy_str
B_w_s = sigy_wk/sigx_str
xi=1#ksi
ho =1# ksi


# %%


def DQx_DQy_L2(bb_name,
            ax_,
            ay_,
            dx_sig,
            dy_sig,
            A_w_s,
            B_w_s,
            r,
            xi,
            ho):
   return {bb_name:  (dtuneL2.DQx_DQy(ax_,
                                      ay_,
                                      dx_sig,
                                      dy_sig,
                                      A_w_s,
                                      B_w_s,
                                      r,
                                      xi,
                                      ho,
                                      Model))}

# %%

axtab=np.array(ax)
aytab=np.array(ay)
result=DQx_DQy_L2(
    bbname,
    axtab,
    aytab,
    dx_sig,
    dy_sig,
    A_w_s,
    B_w_s,
    r,
    xi,
    ho
    )


# %%
# analytic tune shifts
dQx,dQy=result[bbname][0] ,result[bbname][1]

dQx=dQx*ksi
dQy=dQy*ksi

dQx=dQx+Qx_0-np.rint(Qx_0)
dQy=dQy+Qy_0-np.rint(Qy_0)
# %%
plt.plot(dQx,dQy, ls="",    
         marker='x',    
         markersize=10,    
         markerfacecolor="none",    
         markeredgecolor="red",    
         markeredgewidth=2,    
         label="analytic")
plt.grid()
plt.axis('square')
plt.title('Model='+Model+ " ksi="+str(round(ksi,6)))

def f_part(x):
    return modf(x)[1]
xtxt=r'$\Delta Q_x^{BBLR}[\xi]$'
ytxt=r'$\Delta Q_y^{BBLR}[\xi$'
plt.xlabel(xtxt)
plt.ylabel(ytxt)
plt.grid()

do_get_foot=True
if do_get_foot:
#    plt.figure(1)
    fp0 = mini.get_footprint(nemitt_x=nemitt_x, nemitt_y=nemitt_y,
                         freeze_longitudinal=True)
    fp0.plot(color='k', label='get_footprint')
##############

plt.grid()
plt.axis('square')



# %%
