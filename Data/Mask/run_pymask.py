import subprocess
import sys
import os


if 'BBStudies/Executables/py_BB/bin' not in os.environ.get('PATH').split(':')[0]:
    raise Exception('Wrong Distribution')

# Running pymask
cwd = os.getcwd()
os.chdir('./')
exec(open("000_pymask_rich.py").read())
opticsFile = configuration['optics_file'].split('/')[-1]
os.chdir(cwd)

# NOTE: Make sure you add: 
# pm.install_lenses_in_sequence(mad_tra[]ck, bb_dfs['b2'], 'lhcb2')
# at line 438 of '000_pymask_rich.py'{}
# Saving sequences and BB dfs
for seq in ['lhcb1','lhcb2']:
    mad_track.input(f'use, sequence={seq};')
    mad_track.twiss()
    mad_track.survey()
    
    twiss = mad_track.table.twiss.dframe()
    survey = mad_track.table.survey.dframe()

    twiss.to_pickle(f"../Checks/{seq}_opticsfile{opticsFile.split('.')[-1]}_twiss.pkl")
    survey.to_pickle(f"../Checks/{seq}_opticsfile{opticsFile.split('.')[-1]}_survey.pkl")
    
bb_dfs['b1'].to_pickle(f"../Checks/lhcb1_opticsfile{opticsFile.split('.')[-1]}_bb_dfs.pkl")
bb_dfs['b2'].to_pickle(f"../Checks/lhcb2_opticsfile{opticsFile.split('.')[-1]}_bb_dfs.pkl")