
# Choose the list of IRs
IR_list=['1','5']
#IR_list=['1']
#IR_list=['1','5','2','8']
print("Using IR_list=",IR_list)

# Choose the collision Model (Beam-Beam Long-Range, Ideal Wire or OCTupole)  
# Note: Both BBLR and HO collisions are treated with the BBLR Model

Model='BBLR'
#Model='IW'
Model='OCT'

# If one wants to use only HO or Long Range. If not, both False.  
OnlyHo=False
OnlyLR=False

############################
import numpy as np
import sys
import xtrack as xt

sys.path.append('../')

import matplotlib.pyplot as plt

import BBStudies.Physics.Detuning_L2 as dtuneL2
from scipy import constants

import BBStudies.Physics.Base as phys
import BBStudies.Plotting.BBPlots as bbplt
import BBStudies.Physics.Constants as cst
import scipy.stats as sciStat
from multiprocessing import Pool
import time

## %%
#### changed to in in 1_build coll changed optics 
### opticsfile.34 change 
# in confg.yaml 
my_json = ('../BBStudies/Run3_configuration/'
           '2_configure_and_track/'
           'final_collider.json')
collider = xt.Multiline.from_json(my_json)
collider.build_trackers()
## %%


## %%
# Here a_min =.5 so that using ERR=1.e-5 in Detuning_L2 is accurate enough.
# Very small amplitudes min_a << 1 may require better integration accuracy )
a_min = .5
a_max = 8
n_amp = 15
n_wings = 16
coordinates = phys.polar_grid(r_sig=np.linspace(a_min, a_max, n_amp),
                              theta_sig=np.linspace(
                                  0.05*np.pi/2, 0.95*np.pi/2, n_wings),
                              emitt=[2.5e-6/7000, 2.5e-6/7000])
ax = coordinates['x_sig']
ay = coordinates['y_sig']

energy = 7000
# gamma relativistic of a proton at 7 TeV
gamma_rel = energy / \
    (constants.physical_constants['proton mass energy equivalent in MeV'][0]/1000)
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


# Compute dummy xi
Nb = 1e11
round_emitt = 2.5e-6
xi = Nb*cst.r_p/(4*np.pi*round_emitt)
number_of_ho_slices = 11

## %%
def DQx_DQy_L2(bb_name,  r,
                dx_sig,
                dy_sig,
                A_w_s,
                B_w_s,
                xi,
                ho):
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



    ## %%
 
 
    my_filter_string = 'bb_(ho|lr)\.(r|l|c)'+IR+'.*'
    if OnlyHo:
        my_filter_string = 'bb_(ho)\.(r|l|c)'+IR+'.*'
    elif OnlyLR:
        my_filter_string = 'bb_(lr)\.(r|l|c)'+IR+'.*'

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

    fw = 1
    r = sigma_y_strong/sigma_x_strong

    name_weak = twiss_filtered[beam_weak].name


    ## %%
    xi_list = []
    ho_list = []
    for my_index, bb in enumerate(name_weak):
        if bb.find('bb_ho') != -1:
            xi_list.append(xi/number_of_ho_slices)
            ho_list.append(True)
        else:
            xi_list.append(xi)
            ho_list.append(False)
 

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
                                               ho_list))
    e_time = time.time()
    print(f'Execution time for IR{IR}, {(e_time-s_time):.3f} s')

    ## %%

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


Qx_0,Qy_0 = 0.31, 0.32
DQx_All = 0
DQy_All = 0
for _IR in IR_list:
    DQx_IR,DQy_IR=Get_Footprint(IR=_IR) 
    DQx_All += DQx_IR
    DQy_All += DQy_IR
    ##plt.plot(DQx_All , DQy_All ,'xr')
    ##plt.grid()
    #plt.axis('square')

fp_x = DQx_All + Qx_0
fp_y = DQy_All + Qy_0

for window in [0.01,0.05]:

    Qx_lim    = [Qx_0-4*window/4,Qx_0+window/4]
    Qy_lim    = [Qy_0-4*window/4,Qy_0+window/4]

    plt.figure(figsize=(6,6))
    bbplt.workingDiagram(order=12,Qx_range=Qx_lim,Qy_range = Qy_lim,alpha=0.15)

    bbplt.polarmesh(fp_x,fp_y,alpha=0.1,r=coordinates['r_sig'],
                    theta=coordinates['theta_sig'],color='darkslateblue')
    plt.scatter(fp_x,fp_y,s = 30*sciStat.norm.pdf(coordinates['r_sig'])
                /sciStat.norm.pdf(0),zorder=10)
    plt.plot(Qx_0,Qy_0,'P',color='C3',alpha=0.8,label='Unperturbed')

    plt.legend(loc='upper right')
    plt.axis('square')
    plt.xlim(Qx_lim)
    plt.ylim(Qy_lim)
    plt.tight_layout()
plt.show()

plt.plot(ax, ay, 'o')
plt.xlabel('ax [$\sigma$]')
plt.ylabel('ay [$\sigma$]')
plt.axis('square')
plt.close()
# normalized footprint must agree with old Dobrin_Example_L2.py
plt.plot(DQx_All/xi, DQy_All/xi,
         ls="--",marker="x",
         label=IR_list, markersize=4,color='black' )
plt.grid()
plt.axis('square')
plt.title('Model='+Model)





