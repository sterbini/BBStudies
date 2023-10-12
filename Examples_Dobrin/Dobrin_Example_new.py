# %%
import numpy as np
import sys
import xtrack as xt

sys.path.append('../')
import BBStudies.Physics.Detuning as dtune

import BBStudies.Physics.Base as phys
import BBStudies.Plotting.BBPlots as bbplt
import BBStudies.Physics.Constants as cst

# %%
my_json = ('../BBStudies/Run3_configuration/'
           '2_configure_and_track/'
           'final_collider.json')
collider = xt.Multiline.from_json(my_json)
collider.build_trackers()


line   = {}
twiss  = {}
survey = {}

#collider.vars['vrf400'] = 12.0
#collider.vars['lagrf400.b1'] = 0.5
#collider.vars['lagrf400.b2'] = 0.0

# Twiss and survey from xsuite line
for seq in ['lhcb1','lhcb2']:
    _beam = seq[-2:]
    
    # Importing Line
    
    # tracker -> twiss + survey
    if _beam == 'b1':
        twiss[_beam]  = collider[seq].twiss()#.to_pandas(index="name")
        survey[_beam] = collider[seq].survey(element0='ip1')#.to_pandas(index="name")
    elif _beam == 'b2':
        twiss[_beam]   = collider[seq].twiss().reverse()#.to_pandas(index="name")
        survey[_beam]  = collider[seq].survey(element0='ip1').reverse()#.to_pandas(index="name")


assert (survey['b1']['X','ip1']-survey['b2']['X','ip1'])==0

# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sciStat
import sys

sys.path.append('../')
import BBStudies.Physics.Detuning as dtune
import BBStudies.Physics.Base as phys
import BBStudies.Physics.Constants as cst
import BBStudies.Plotting.BBPlots as bbplt


coordinates = phys.polar_grid(  r_sig     = np.linspace(1,6.5,5),
                                theta_sig = np.linspace(0.05*np.pi/2,0.95*np.pi/2,5),
                                emitt     = [2.5e-6/7000,2.5e-6/7000])


# %%
twiss_filtered = {}
survey_filtered = {}

my_filter_string = 'bb_(ho|lr)\.(r|l|c)1.*'
#Xmy_filter_string = 'bb_(ho)\.(r|l|c)1.*'

#my_filter_string = 'bb_(lr)\.(r|l|c)(1|5).*25'

for beam in ['b1','b2']:
    twiss_filtered[beam]  = twiss[beam][:, my_filter_string]
    survey_filtered[beam]  = survey[beam][['X','Y','Z'], my_filter_string]


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
    twiss_filtered[beam_weak]['x'] - twiss_filtered[beam_strong]['x'] +
    survey_filtered[beam_weak]['X']- survey_filtered[beam_strong]['X']
    )
d_y_weak_strong_in_meter = (
    twiss_filtered[beam_weak]['y'] - twiss_filtered[beam_strong]['y'] +
    survey_filtered[beam_weak]['Y']- survey_filtered[beam_strong]['Y']
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

# %%
plt.plot(ax,ay,'o')
plt.xlabel('ax [$\sigma$]')
plt.ylabel('ay [$\sigma$]')
plt.axis('square')

# %%

plt.plot(s,np.abs(dx_sig),'o')
plt.xlabel('s [m]') 
plt.ylabel('distance in x [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)
# %%
plt.plot(s, np.abs(dy_sig),'o')
plt.xlabel('s [m]') 
plt.ylabel('distance in y [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)
## %%
#plt.plot(s, r,'o')
#plt.xlabel('s [m]') 
#plt.ylabel('r')
#plt.title('Filtering by '+my_filter_string)
# %%
#plt.plot(s, A_w_s,'o')
#plt.xlabel('s [m]') 
#plt.ylabel('A_w_s')
#plt.title('Filtering by '+my_filter_string)
# %%
#plt.plot(s, B_w_s,'o')
#plt.xlabel('s [m]') 
#plt.ylabel('B_w_s')
#plt.title('Filtering by '+my_filter_string)

# %%
""""
DQx_HO,DQy_HO = dtune.DQx_DQy(ax, 
              ay, 0, 
                      0, 
                      1,  
                      1,  
                      1,  
                      xi=1e-2, 
                      fw=1)
"""
# %%
# import time
# import BBStudies.Physics.Constants as cst

# # Initialize tuneshift
# DQx_HO,DQy_HO = np.zeros(len(ax)),np.zeros(len(ax))



# #  %%

# # Sum tuneshift for all head-on
# for my_index, bb in enumerate(name_weak):
#     print(f'Computing {bb}...')
#     print(f'Index: {my_index}')

#     s_time = time.time()

#     _DQx,_DQy = dtune.DQx_DQy(  ax     = coordinates['x_sig'],
#                                 ay     = coordinates['y_sig'],
#                                 r      = r[my_index],
#                                 dx_sig = dx_sig[my_index],
#                                 dy_sig = dy_sig[my_index],
#                                 A_w_s  = A_w_s[my_index],
#                                 B_w_s  = B_w_s[my_index],
#                                 xi     = xi_list[my_index])

#     DQx_HO += _DQx
#     DQy_HO += _DQy

#     e_time = time.time()
#     print(f'Execution time, {bb}: {(e_time-s_time):.3f} s')


# %%

# # Close-up and zoomed out plot
# #====================================================
Qx_0,Qy_0 = 0.31, 0.32
#fp_x = Qx_0 + DQx_HO
#fp_y = Qy_0 + DQy_HO
#plt.plot(fp_x,fp_y,'o') 

# #====================================================
 
# %%
from multiprocessing import Pool

# Compute dummy xi
Nb = 1e11
round_emitt = 2.5e-6
xi = Nb*cst.r_p/(4*np.pi*round_emitt)
number_of_ho_slices = 11


xi_list = []
for my_index, bb in enumerate(name_weak):
    if bb.find('bb_ho') != -1:
        xi_list.append(xi/number_of_ho_slices)
    else:
        xi_list.append(xi)
def func(a, b):
    return {'test':(a + b, a * b)}



def myDQx_DQy(bb_name,  r,
                dx_sig, 
                dy_sig, 
                A_w_s, 
                B_w_s, 
                xi ):
    return {bb_name: (dtune.DQx_DQy(coordinates['x_sig'],
                                    coordinates['y_sig'],
                                dx_sig,
                                dy_sig,
                                A_w_s ,
                                B_w_s ,
                                r     ,
                                xi,
                                fw=1   ))}

#myDQx_DQy(name_weak[0],r[0],dx_sig[0],dy_sig[0],A_w_s[0],B_w_s[0],xi_list[0])
it=3

result1=myDQx_DQy(name_weak[it],r[it],dx_sig[it],dy_sig[it],A_w_s[it],B_w_s[it],xi_list[it])
key_name,=result1
# %%
import BBStudies.Physics.Detuning_L2 as dtuneL2
#UseModel='IW'
def DQx_DQy_dob(bb_name,  r,
                dx_sig, 
                dy_sig, 
                A_w_s, 
                B_w_s,
                xi):
   return {bb_name:  (dtuneL2.DQx_DQy(coordinates['x_sig'],
                                     coordinates['y_sig'],
                                dx_sig,
                                dy_sig,
                                A_w_s ,
                                B_w_s ,
                                r     ,
                                xi ))}    
    
result_no_nan=DQx_DQy_dob(name_weak[it],r[it],dx_sig[it],dy_sig[it],A_w_s[it],B_w_s[it],xi_list[it])

 
# %%
print(" compare the two functions DQx_DQy -- old and new \n is this zero? \n ",(result1[key_name][0]-result_no_nan[key_name][0]))


# %%
with Pool(64) as pool:
    result = pool.starmap(DQx_DQy_dob, zip(name_weak, 
                                        r, 
                                        dx_sig, 
                                        dy_sig, 
                                        A_w_s, 
                                        B_w_s, 
                                        xi_list))

# %%  # former BUG
"""
dtune.DQx_DQy(  ax     = [coordinates['x_sig'].values[1]],
                                ay     = [coordinates['y_sig'].values[1]],
                                r      = r[0],
                                dx_sig = dx_sig[0],
                                dy_sig = dy_sig[0],
                                A_w_s  = A_w_s[0],
                                B_w_s  = B_w_s[0],
                                xi     = xi_list[0])


"""
# %%

# convert a list of dict in a dict
dict_result = {}

for my_dict in result:
    for key in my_dict:
        dict_result[key] = my_dict[key]

# define a function that takes a np.array and replace nan with 0
def nan_to_zero(my_array):
    my_array[np.isnan(my_array)] = 0
    return my_array

DQx_HO= 0
DQy_HO = 0
for my_bb in dict_result.keys():
    print(my_bb)
    DQx_HO+=nan_to_zero(dict_result[my_bb][0])
    DQy_HO+=nan_to_zero(dict_result[my_bb][1])

# %%



Qx_0,Qy_0 = 0.31, 0.32

fp_x = Qx_0 + DQx_HO
fp_y = Qy_0 + DQy_HO
 
# %%
plt.plot(fp_x,fp_y,'--xr')

plt.axis('square')


# %%
