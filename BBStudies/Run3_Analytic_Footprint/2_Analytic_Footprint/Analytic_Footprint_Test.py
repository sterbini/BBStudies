# %%
"""
The collider built in 1_build_collider (round or flat optics) is first tracked for
footprint, then the analytic footprint is computed and compared.
Dobrin, Dec 2023
"""
# ==================================================================================================
# --- Imports
# ==================================================================================================
import json
import yaml
import time
import logging
import numpy as np
from numpy import *
import pandas as pd
import os

import matplotlib.pyplot as plt
import numpy as np
from numpy import *

from scipy import constants

from scipy.optimize import minimize_scalar
import xtrack as xt
import tree_maker

import sys
#sys.path.append('/home1/BBStudies')

sys.path.append('../') 
sys.path.append('../../../')

#import pprint
#pprint.pprint( sys.path)

from misc import generate_orbit_correction_setup
from misc import luminosity_leveling, luminosity_leveling_ip1_5


#run3_tools =__import__('BBStudies.Run3_configuration.2_configure_and_track')

import BBStudies.Physics.Base as phys

from copy_of_2_configure_and_track import *

import BBStudies.Physics.Detuning_L2  as dtuneL2

#import os
#print(os.getcwd())

# %%

# ==================================================================================================
# --- Define local function for collider configuration (flag to switch matching off)
# ==================================================================================================
def configure_collider_no_match(
    config_sim,
    config_collider,
    skip_beam_beam=False,
    match_tune=False,
    save_collider=False,
    return_collider_before_bb=False,
):
    # Generate configuration files for orbit correction
    generate_configuration_correction_files()

    # Rebuild collider
    collider = xt.Multiline.from_json(config_sim["collider_file"])

    # Install beam-beam
    collider, config_bb = install_beam_beam(collider, config_collider)

    # Build trackers
    collider.build_trackers()

    # Set knobs
    collider, conf_knobs_and_tuning = set_knobs(config_collider, collider)

    # Match tune and chromaticity
    if match_tune:
        collider = match_tune_and_chroma(
            collider, conf_knobs_and_tuning, match_linear_coupling_to_zero=False
        )
        
    # Compute the number of collisions in the different IPs
    n_collisions_ip1_and_5, n_collisions_ip2, n_collisions_ip8 = compute_collision_from_scheme(
        config_bb
    )

    # Do the leveling if requested
    if "config_lumi_leveling" in config_collider and not config_collider["skip_leveling"]:
        collider = do_levelling(
            config_collider, config_bb, n_collisions_ip8, collider, n_collisions_ip1_and_5
        )
    else:
        print(
            "No leveling is done as no configuration has been provided, or skip_leveling"
            " is set to True."
        )

    # Add linear coupling
    #collider = add_linear_coupling(conf_knobs_and_tuning, collider)

    # Rematch tune and chromaticity
    if match_tune:
        collider = match_tune_and_chroma(
        collider, conf_knobs_and_tuning, match_linear_coupling_to_zero=True
    )
#
#    # Assert that tune, chromaticity and linear coupling are correct one last time
#    assert_tune_chroma_coupling(collider, conf_knobs_and_tuning)

    # Return twiss and survey before beam-beam if requested
    if return_collider_before_bb:
        print("Saving collider before beam-beam configuration")
        collider_before_bb = xt.Multiline.from_dict(collider.to_dict())

    if not skip_beam_beam:
        # Configure beam-beam
        collider = configure_beam_beam(collider, config_bb)

    if save_collider:
        # Save the final collider before tracking
        collider.to_json("final_collider.json")

    if return_collider_before_bb:
        return collider, config_bb, collider_before_bb
    else:
        return collider, config_bb

# ==================================================================================================
# --- Function to set all high orders to zero
# ==================================================================================================

def Set_High_Ord_Mult_Zero(collider):
    for ee, nn in zip(collider['lhcb1'].elements, collider['lhcb1'].element_names):
    #    print(ee.__class__.__name__ =="Multipole" )
        if ee.__class__.__name__ =="Multipole":
            kk=ee.knl
            if len(kk)>2:
                mult=kk[1:]
                if linalg.norm(mult)>0:
                    ee.knl=0
    # set all skew to zero
    for ee, nn in zip(collider['lhcb1'].elements, collider['lhcb1'].element_names):
    #    print(ee.__class__.__name__ =="Multipole" )
        if ee.__class__.__name__ =="Multipole":
            kk=ee.ksl
            if len(kk)>2:
                mult=kk[1:]
                if linalg.norm(mult)>0:
                    ee.ksl=0
    
    
    for ee, nn in zip(collider['lhcb1'].elements, collider['lhcb1'].element_names):
        if ee.__class__.__name__ =="Multipole":
            kk=ee.knl
            #kk=ee.ksl
            if len(kk)>2:
                if np.linalg.norm(kk)>0:
                    print("not zero yet ",kk)
    



# %%
##########################################################################################
# ...        Collider with and without beam-beam 
#            Choose round ot flat optics from ../1_build_collider ..    
#########################################################################################

config_path="config.yaml"
config, config_sim, config_collider = read_configuration(config_path)

use_flat_opt=True #False
if use_flat_opt:
    config_sim["collider_file"]='../1_build_collider/collider_flat/collider.json'
else:
    config_sim["collider_file"]='../1_build_collider/collider/collider.json'


collider, config_bb,collider_before_bb= configure_collider_no_match(
config_sim, 
config_collider, 
skip_beam_beam=False,
match_tune=False,
return_collider_before_bb=True,
save_collider=False
)

collider_before_bb.build_trackers()

Set_High_Ord_Mult_Zero(collider_before_bb)

Set_High_Ord_Mult_Zero(collider)

# %%
# =================================================================================
# ---  Collider tunes before BB: Qx0,Qy0
# =================================================================================

tw_bef_bb = collider_before_bb['lhcb1'].twiss(method='4d')

# tune and chromaticity values
print(f"Qx = {tw_bef_bb['qx']:.5f} Qy = {tw_bef_bb['qy']:.5f} "
      f"Q'x = {tw_bef_bb['dqx']:.5f} Q'y = {tw_bef_bb['dqy']:.5f}")


print("beta star IP1 no bb")
flat_opt_check=tw_bef_bb[:,'ip1'].betx[0]/tw_bef_bb[:,'ip1'].bety[0]
print(tw_bef_bb[:,'ip1'].betx[0])
print(tw_bef_bb[:,'ip1'].bety[0])


Qx_0=tw_bef_bb['qx']    
Qy_0=tw_bef_bb['qy']    

Qx0,Qy0=Qx_0-np.rint(Qx_0),Qy_0-np.rint(Qy_0)


# %%
# =================================================================================
# --- Use same emittance in both planes (but possibly flat optics)
# =================================================================================
# 
energy=7000

nemit=config_collider['config_beambeam']['nemitt_x']

Nb=config_collider['config_beambeam']['num_particles_per_bunch']

number_of_ho_slices =config_collider['config_beambeam']['num_slices_head_on']

emittance_strong_nx=nemit
emittance_strong_ny=nemit
emittance_weak_nx=nemit
emittance_weak_ny=nemit
gamma_rel = energy / (constants.physical_constants['proton mass energy equivalent in MeV'][0]/1000)
beta_rel = np.sqrt(1-1/gamma_rel**2)
emittance_strong_x = emittance_strong_nx/gamma_rel/beta_rel
emittance_strong_y = emittance_strong_ny/gamma_rel/beta_rel
emittance_weak_x = emittance_weak_nx/gamma_rel/beta_rel
emittance_weak_y = emittance_weak_ny/gamma_rel/beta_rel

rp = 1.534698E-18
ksi= (Nb*rp)/(4*np.pi*nemit)
ksi_txt=str(round(ksi,5))

print('Beam-beam parameter xi=',ksi_txt)

# %%
# =================================================================================
# --- Define radial uniform amplitude distribution. Use the standard get_footprint. 
# =================================================================================
# 
a_min = .05
a_max = 8
n_amp = 8
n_wings =8
n_turns=1024


fp0 = collider['lhcb1'].get_footprint(nemitt_x=nemit, nemitt_y=nemit,
                 freeze_longitudinal=True,r_range=(a_min,a_max),
                 n_r=n_amp,
                 n_theta=n_wings,
                 n_turns=n_turns )
#fp0.plot(color='k', label='get_footprint')

coordinates = phys.polar_grid(r_sig=np.linspace(a_min, a_max, n_amp),
                              theta_sig=np.linspace(
                                  0.05*np.pi/2, 0.95*np.pi/2, n_wings),
                              emitt=[emittance_weak_x,emittance_weak_y])

ax = coordinates['x_sig']
ay = coordinates['y_sig']
ax=np.array(ax)
ay=np.array(ay)

fig=plt.figure(figsize=(6, 3))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
fpx=ndarray.flatten(fp0.qx) 
fpy=ndarray.flatten(fp0.qy) 
ax1.scatter(ax,ay,c='black',ls='')

ax2.plot(fpx,fpy, ls="",    
         marker='o',    
         markersize=10,    
         markerfacecolor="none",    
         markeredgecolor="green",    
         markeredgewidth=1,    
         label="get foot")
ax1.axis('square')
ax2.axis('square')
ax1.grid()
ax2.grid()

plt.suptitle("get_footprint, green "+" ksi="+ksi_txt)
plt.tight_layout()



# %%
# ==========================================================================
# --- Analytic calculation of footrpint with 2D-Bessel - script Detuning_L2 
#     Choose the collision Model (Beam-Beam Long-Range, Ideal Wire or OCTupole)  
#     Note:  BBLR and HO collisions are treated with BBLR 
#================================================================================
Model = 'BBLR'
#Model = 'IW'
#Model = 'OCT'
from multiprocessing import Pool
import time

IR_list=['1','5','2','8']
def DQx_DQy_L2(bb_name,  r,
                dx_sig,
                dy_sig,
                A_w_s,
                B_w_s,
                xi,
                ho,
                state):
    if state == 0:
        return {bb_name: [0,0]}
    else:
        return {bb_name:  (dtuneL2.DQx_DQy(coordinates['x_sig'],
                                      coordinates['y_sig'],
                                      dx_sig,
                                      dy_sig,
                                      A_w_s,
                                      B_w_s,
                                      r,
                                      xi,
                                      ho,
                                      Model))}


def Get_Footprint(IR):
    print('computing IR'+IR)
    twiss = {}
    survey = {}
    # Twiss and survey from xsuite line
    for seq in ['lhcb1', 'lhcb2']:
        _beam = seq[-2:]

        # Importing Line
        # tracker -> twiss + survey
        if _beam == 'b1':
            twiss[_beam] = collider[seq].twiss()   
            survey[_beam] = collider[seq].survey(
                element0='ip'+IR)   
        elif _beam == 'b2':             
            twiss[_beam] = collider[seq].twiss().reverse()
            survey[_beam] = collider[seq].survey(
                element0='ip'+IR).reverse()   

    my_filter_string = 'bb_(ho|lr)\.(r|l|c)'+IR+'.*'
    twiss_filtered = {}
    survey_filtered = {}
    for beam in ['b1', 'b2']:
        twiss_filtered[beam] = twiss[beam][:, my_filter_string]
        survey_filtered[beam] = survey[beam][['X', 'Y', 'Z'], my_filter_string]


    beam_weak = 'b1'
    beam_strong = 'b2'

    s = survey_filtered[beam_strong]['Z']
    d_x_weak_strong_in_meter = (
        twiss_filtered[beam_weak]['x'] - twiss_filtered[beam_strong]['x'] +
        survey_filtered[beam_weak]['X'] - survey_filtered[beam_strong]['X']
    )
    d_y_weak_strong_in_meter = (
        twiss_filtered[beam_weak]['y'] - twiss_filtered[beam_strong]['y'] +
        survey_filtered[beam_weak]['Y'] - survey_filtered[beam_strong]['Y']
    )

    sigma_x_strong = np.sqrt(
        twiss_filtered[beam_strong]['betx']*emittance_strong_x)
    sigma_y_strong = np.sqrt(
        twiss_filtered[beam_strong]['bety']*emittance_strong_y)

    sigma_x_weak = np.sqrt(twiss_filtered[beam_weak]['betx']*emittance_weak_x)
    sigma_y_weak = np.sqrt(twiss_filtered[beam_weak]['bety']*emittance_weak_y)

    dx_sig = d_x_weak_strong_in_meter/sigma_x_strong
    dy_sig = d_y_weak_strong_in_meter/sigma_y_strong

    A_w_s = sigma_x_weak/sigma_y_strong
    B_w_s = sigma_y_weak/sigma_x_strong
    r = sigma_y_strong/sigma_x_strong

    name_weak = twiss_filtered[beam_weak].name

    xi_list = []
    ho_list = []
    count_ho=0
    for my_index, bb in enumerate(name_weak):
        if bb.find('bb_ho') != -1:
            xi_list.append(ksi/number_of_ho_slices)
            ho_list.append(True)
            count_ho=count_ho+1
        else:
            xi_list.append(ksi)
            ho_list.append(False)
    assert config_bb["num_slices_head_on"] == count_ho
    #state true means active bunch 
    state_list=[]
    for my_index, bb in enumerate(name_weak):
        state=collider.vars[bb + '_scale_strength']._value
        state_list.append(state)

    ## %%
    s_time = time.time()
    with Pool(64) as pool:
        result = pool.starmap(DQx_DQy_L2, zip(name_weak,
                                               r,
                                               dx_sig,
                                               dy_sig,
                                               A_w_s,
                                               B_w_s,
                                               xi_list,
                                               ho_list,
                                               state_list))
    e_time = time.time()
    print(len(name_weak))

    assert 2*config_bb["num_long_range_encounters_per_side"]['ip'+IR] == (len(name_weak)-count_ho)

    print(f'Execution time for IR{IR}, {(e_time-s_time):.3f} s')

    # convert a list of dict in a dict
    dict_result = {}

    for my_dict in result:
        for key in my_dict:
            dict_result[key] = my_dict[key]

    ## %%

    DQx_IR = 0
    DQy_IR = 0
    for my_bb in dict_result.keys():
    #    print(my_bb)
        DQx_IR += dict_result[my_bb][0]
        DQy_IR += dict_result[my_bb][1]

    return DQx_IR, DQy_IR
# %%
DQx_All = 0
DQy_All = 0
for _IR in IR_list:
    DQx_IR,DQy_IR=Get_Footprint(IR=_IR) 
    plt.plot(DQx_IR/ksi,DQy_IR/ksi, ls="-",   
         marker='x',    
         label="IR"+_IR)
    DQx_All += DQx_IR
    DQy_All += DQy_IR
plt.grid()
plt.axis('square')
plt.legend()


# %%
####################################################################
# --- Plot analytic and numerical footprints
#####################################################################

show_in_xi_units=True #False

dQx=DQx_All   
dQy=DQy_All  
if show_in_xi_units:
    dQx=dQx/ksi
    dQy=dQy/ksi
plt.plot(dQx,dQy, ls="",    
         marker='x',    
         markersize=10,    
         markerfacecolor="none",    
         markeredgecolor="red",    
         markeredgewidth=2,    
         label="analytic")
plt.grid()
plt.axis('square')

fpx=ndarray.flatten(fp0.qx)-Qx0 
fpy=ndarray.flatten(fp0.qy)-Qy0

if show_in_xi_units:
    fpx=fpx/ksi
    fpy=fpy/ksi

plt.plot(fpx,fpy, ls="",    
         marker='o',    
         markersize=10,    
         markerfacecolor="none",    
         markeredgecolor="green",    
         markeredgewidth=2,    
         label="get_footprint")

                
plt.grid()
plt.axis('square')

if show_in_xi_units:
    xtxt=r'$\Delta Q_x^{BBLR}[\xi]$'
    ytxt=r'$\Delta Q_y^{BBLR}[\xi]$'
else:
    xtxt=r'$\Delta Q_x^{BBLR}$'
    ytxt=r'$\Delta Q_y^{BBLR}$'
     
plt.xlabel(xtxt)
plt.ylabel(ytxt)
plt.grid()
plt.title('Model='+Model+ 
r'  $\xi = $' + str(round(ksi,5))+" i_bunch_b1 ="+
str(config_collider['config_beambeam']['mask_with_filling_pattern']['i_bunch_b1'])+ "\n"+
str(config_collider['config_beambeam']['mask_with_filling_pattern']['pattern_fname'] )+
"\n"+
r'  $a_{min} = $' + str(round(a_min,2))+
r'  $a_{max} = $' + str(round(a_max,2))+
           "\n"
r'$\beta_x^{IP1} / \beta_y^{IP1} =$'+str(round(flat_opt_check,2))
 )
plt.legend()
 
# %%
