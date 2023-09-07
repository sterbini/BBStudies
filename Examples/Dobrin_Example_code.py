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
        survey[_beam]  = collider[seq].survey(element0='ip1').reverse() #.to_pandas(index="name")


# # Twiss and survey from xsuite line
# for seq in ['lhcb1','lhcb2']:
#     _beam = seq[-2:]
    
#     # Importing Line
    
#     # tracker -> twiss + survey
#     if _beam == 'b1':
#         twiss[_beam]  = collider[seq].twiss().reverse()#.to_pandas(index="name")
#         survey[_beam] = collider[seq].survey(element0='ip1').reverse()#.to_pandas(index="name")
#     elif _beam == 'b2':
#         twiss[_beam]   = collider[seq].twiss()#.to_pandas(index="name")
#         survey[_beam]  = collider[seq].survey(element0='ip1') #.to_pandas(index="name")



# %%
#print(survey['b1']['X','ip1'],survey['b2']['X','ip1'])
print( "IP1 is HOR  xing ",twiss['b1']['px','ip1'])
print( "IP5 is VERT xing ",twiss['b1']['py','ip5'])

# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sciStat
import sys

#sys.path.append('../')
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
my_filter_string = 'bb_(ho|lr)\.(r|l|c)1.*'
my_filter_string = 'bb_(lr)\.(r|l|c)1.*'
#dobrin choose only slot #12
my_filter_string = 'bb_(lr)\.(r|l|c)(1|5).*_12'

#dobrin choose only slot #12 5
my_filter_string = 'bb_(lr)\.(r|l|c)1.*12'
my_filter_string = 'bb_(lr)\.(r|l|c)1.*'

#dobrin choose only slot outside #.. 
#my_filter_string = 'bb_(lr)\.(r|l|c)1.*22'

for beam in ['b1','b2']:
    twiss_filtered[beam]  = twiss[beam][:, my_filter_string]
    survey_filtered[beam]  = survey[beam][['X','Y','Z'], my_filter_string]
names=twiss_filtered['b1']['name']
 
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
print(d_x_weak_strong_in_meter, d_y_weak_strong_in_meter)
print(dx_sig, dy_sig,r)  

# %%

plt.plot(s,(dx_sig),'o')
plt.xlabel('s [m]') 
plt.ylabel('sep in x [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)

# %%
plt.plot(s, (dy_sig),'o')
plt.xlabel('s [m]') 
plt.ylabel('sep in y [$\sigma$ of the strong beam]')
plt.title('Filtering by '+my_filter_string)

# %%
plt.plot(s, r,'o')
plt.xlabel('s [m]') 
plt.ylabel('r')
plt.title('Filtering by '+my_filter_string)
# %%
plt.plot(s, A_w_s,'o')
plt.xlabel('s [m]') 
plt.ylabel('A_w_s')
plt.title('Filtering by '+my_filter_string)

# %%
plt.plot(s, B_w_s,'o')
plt.xlabel('s [m]') 
plt.ylabel('B_w_s')
plt.title('Filtering by '+my_filter_string)

# %%

# %%

import time
import BBStudies.Physics.Constants as cst

# Initialize tuneshift
DQx_tot,DQy_tot = np.zeros(len(coordinates)),np.zeros(len(coordinates))

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
print("bb parameters are xi=",xi_list)        

model=0
# Sum tuneshift for all
for my_index, bb in enumerate(name_weak):
    print(f'Computing {bb}...')
    print(f'Index: {my_index}')

    s_time = time.time()
    _dx_sig = dx_sig[my_index]
    _dy_sig = dy_sig[my_index]
    sep_reset=20
    print(_dx_sig,_dy_sig)
    
    """
    if np.abs(_dx_sig)<0.1 :               
        _dx_sig=0.0
    if np.abs(_dy_sig)<0.1 :
        _dy_sig=0.0
        
    
    if np.abs(_dx_sig)>sep_reset :
        print('dx=',_dx_sig )
        print("reset dx to "+ str(sep_reset))
        _dx_sig=sep_reset*np.sign(_dx_sig)
    if np.abs(_dy_sig)>sep_reset :
        print('dy=',_dy_sig)
        print("reset dy to "+str(sep_reset))
        _dy_sig=sep_reset*np.sign(_dy_sig)
    """
        
    _DQx,_DQy = dtune.DQx_DQy_code(  
                                ax     = coordinates['x_sig'],
                                ay     = coordinates['y_sig'],
                                dx_sig =_dx_sig,# dx_sig[my_index],
                                dy_sig =_dy_sig,# dy_sig[my_index],
                                A_w_s  = A_w_s[my_index],
                                B_w_s  = B_w_s[my_index],
                                r      = r[my_index],
                                xi     = xi_list[my_index],
                                 
                                model=model
                                )
   
    DQx_tot += _DQx
    DQy_tot += _DQy

    e_time = time.time()
    print(f'Execution time, {bb}: {(e_time-s_time):.3f} s')

# %%

plt.plot(DQx_tot/xi_list[0],DQy_tot/xi_list[0],
    ls="-",
    marker="x",
    markersize=10,
    markerfacecolor="none",
    markeredgecolor="red",
    markeredgewidth=2,
    label="analytic")
plt.title('model='+str(model)+ ' no reset')
plt.xlabel('DQx [$ xi $]')
plt.ylabel('DQy [$ xi $]')
plt.grid()
plt.axis('square')
 
# %%

# Close-up and zoomed out plot
#====================================================
Qx_0,Qy_0 = 0.31, 0.32
fp_x = Qx_0 + DQx_tot
fp_y = Qy_0 + DQy_tot

for window in [0.02,0.05]:

    Qx_lim    = [Qx_0-4*window/4,Qx_0+window/2]
    Qy_lim    = [Qy_0-4*window/4,Qy_0+window/2]

    plt.figure(figsize=(6,6))
    bbplt.workingDiagram(order=12,Qx_range=Qx_lim,Qy_range = Qy_lim,alpha=0.15)

    bbplt.polarmesh(fp_x,fp_y,alpha=1,r=coordinates['r_sig'],theta=coordinates['theta_sig'],color='darkslateblue')
    plt.scatter(fp_x,fp_y,s = 30*sciStat.norm.pdf(coordinates['r_sig'])/sciStat.norm.pdf(0),zorder=10)
    plt.plot(Qx_0,Qy_0,'P',color='C3',alpha=0.8,label='Unperturbed')

    plt.legend(loc='upper right')
    plt.axis('square')
    plt.xlim(Qx_lim)
    plt.ylim(Qy_lim)
    plt.tight_layout()
#====
# %%
