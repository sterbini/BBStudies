# %%
import numpy as np
import sys
sys.path.append('../')

import xtrack as xt

import BBStudies.Physics.Detuning_code as dtune
import BBStudies.Physics.Base as phys
import BBStudies.Plotting.BBPlots as bbplt
import BBStudies.Physics.Constants as cst


# %%
my_json = ('../BBStudies/Run3_configuration/'
           '2_configure_and_track/'
           'final_collider.json')
collider = xt.Multiline.from_json(my_json)
collider.build_trackers()

# %%
line   = {}
twiss  = {}
survey = {}

# choose line to be saved in mydata -- need IR1 and IR5
ipcase='ip1'
#ipcase='ip5'

#collider.vars['vrf400'] = 12.0
#collider.vars['lagrf400.b1'] = 0.5
#collider.vars['lagrf400.b2'] = 0.0

# Twiss and survey from xsuite line
for seq in ['lhcb1','lhcb2']:
    _beam = seq[-2:]
    
    # Importing Line
     
    # tracker -> twiss + survey
    #  notice that survey is defined wrt element0
    
    if _beam == 'b1':
        twiss[_beam]  = collider[seq].twiss()#.to_pandas(index="name")
        survey[_beam] = collider[seq].survey(element0=ipcase)#.to_pandas(index="name")
    elif _beam == 'b2':
        twiss[_beam]   = collider[seq].twiss().reverse()#.to_pandas(index="name")
        survey[_beam]  = collider[seq].survey(element0=ipcase).reverse() #.to_pandas(index="name")


# %%
print(survey['b1']['X','ip1'],survey['b2']['X','ip1'])
print( "IP1 is HOR  xing ",twiss['b1']['px','ip1'])
print( "IP5 is VERT xing ",twiss['b1']['py','ip5'])

# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sciStat
import sys

sys.path.append('../')
import BBStudies.Physics.Detuning_code as dtune
import BBStudies.Physics.Base as phys
import BBStudies.Physics.Constants as cst
import BBStudies.Plotting.BBPlots as bbplt


coordinates = phys.polar_grid(  r_sig     = np.linspace(1,6.5,5),
                                theta_sig = np.linspace(0.05*np.pi/2,0.95*np.pi/2,6),
                                emitt     = [2.5e-6/7000,2.5e-6/7000])

plt.figure()
plt.plot(coordinates.x_sig,coordinates.y_sig,'o')
plt.axis('square')
# %%
# plt.plot(ax,ay,'o')
# plt.xlabel('ax [$\sigma$]')
# plt.ylabel('ay [$\sigma$]')
# plt.axis('square')

# %%
twiss_filtered = {}
survey_filtered = {}
#my_filter_string = 'bb_(ho|lr)\.(r|l|c)1.*'
#my_filter_string = 'bb_(lr)\.(r|l|c)1.*'
#my_filter_string = 'bb_(lr)\.(r|l|c)1.*'

#dobrin choose only slot #12
#my_filter_string = 'bb_(lr)\.(r|l|c)(1|5).*_12'

#dobrin choose only 
#my_filter_string = 'bb_(lr)\.(r|l|c)1.*'


if ipcase == 'ip1':
	my_filter_string = 'bb_(lr)\.(r|l|c)1.*'	
elif ipcase == 'ip5':
	my_filter_string = 'bb_(lr)\.(r|l|c)5.*'

#dobrin choose only slot outside #.. 
#my_filter_string = 'bb_(lr)\.(r|l|c)1.*12'
#dobrin all 100
#my_filter_string = 'bb_(lr)\.(r|l|c)(1|5).*'


for beam in ['b1','b2']:
    twiss_filtered[beam]  = twiss[beam][:, my_filter_string]
    survey_filtered[beam]  = survey[beam][['X','Y','Z'], my_filter_string]
names=twiss_filtered['b1']['name']
"""
n=len(list(names1))
names=names1
for i in range(n):
    names[i]=names1[i][5:8]+names1[i][10:]
    names[i]=names1[i].upper()
"""
  
 
# %%
from scipy import constants
beam_weak = 'b1'
beam_strong = 'b2'

energy = 7000 
# gamma relativistic of a proton at 7 TeV
gamma_rel = energy/(constants.physical_constants['proton mass energy equivalent in MeV'][0]/1000)
# beta relativistic of a proton at 7 TeV
beta_rel = np.sqrt(1-1/gamma_rel**2)

emittance_strong_nx = 2.5e-6
emittance_strong_ny = 2.5e-6

emittance_weak_nx = 2.5e-6
emittance_weak_ny = 2.5e-6

emittance_strong_x = emittance_strong_nx/gamma_rel/beta_rel
emittance_strong_y = emittance_strong_ny/gamma_rel/beta_rel

emittance_weak_x = emittance_weak_nx/gamma_rel/beta_rel
emittance_weak_y = emittance_weak_ny/gamma_rel/beta_rel

ax = coordinates['x_sig']
ay = coordinates['y_sig']

s = survey_filtered[beam_strong]['Z']
d_x_weak_strong_in_meter = (
    twiss_filtered[beam_weak]['x'] - twiss_filtered[beam_strong]['x'] 
    + survey_filtered[beam_weak]['X']-  survey_filtered[beam_strong]['X']
    )
d_y_weak_strong_in_meter = (
    twiss_filtered[beam_weak]['y'] - twiss_filtered[beam_strong]['y'] 
    +  survey_filtered[beam_weak]['Y']-  survey_filtered[beam_strong]['Y']
    )

sigma_x_strong = np.sqrt(twiss_filtered[beam_strong]['betx']*emittance_strong_x)
sigma_y_strong = np.sqrt(twiss_filtered[beam_strong]['bety']*emittance_strong_y)

sigma_x_weak = np.sqrt(twiss_filtered[beam_weak]['betx']*emittance_weak_x)
sigma_y_weak = np.sqrt(twiss_filtered[beam_weak]['bety']*emittance_weak_y)

dx_sig = d_x_weak_strong_in_meter/sigma_x_strong
dy_sig = d_y_weak_strong_in_meter/sigma_y_strong

A_w_s = sigma_x_weak/sigma_y_strong
B_w_s = sigma_y_weak/sigma_x_strong

fw = 1 
r = sigma_y_strong/sigma_x_strong

name_weak = twiss_filtered[beam_weak].name
#print(d_x_weak_strong_in_meter, d_y_weak_strong_in_meter)
#print(dx_sig, dy_sig,r)  
print(name_weak)

# %%

plt.plot(s,(dx_sig),'o')
plt.ylim(-20, 20)

plt.xlabel('s [m]') 
plt.ylabel('sep in x [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)
plt.xticks(fontsize=11,fontweight = 'bold',rotation = 90)

# %%

plt.plot(s,(dx_sig),'o')
plt.ylim(-20, 20)

plt.xlabel('s [m]') 
plt.ylabel('sep in x [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)
plt.xticks(fontsize=11,fontweight = 'bold',rotation = 90)

# %%

plt.plot(s,(dy_sig),'o')
plt.ylim(-20, 20)

plt.xlabel('s [m]') 
plt.ylabel('sep in y [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)
plt.xticks(fontsize=11,fontweight = 'bold',rotation = 90)

# %%
dir='../Examples_Dobrin/mydata/'

np.savetxt(dir+'names4py_'+ipcase+".dat",np.array(names), fmt='%s')

tmp=np.array([s,dx_sig,dy_sig,r,A_w_s,B_w_s]) 
tmp=np.transpose(tmp)
np.savetxt(dir+'lrtab4py_'+ipcase+".dat",tmp, fmt='%14.10f')

# %%
