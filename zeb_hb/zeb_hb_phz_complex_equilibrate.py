#!/usr/bin/env python
import sys
import os
import shutil

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from gromacs_py import gmx
import pickle

DATA_OUT = 'zeb_hb_phz_complex_sim'

# System Setup
sys_top_folder = os.path.join(DATA_OUT, 'sys_top')
# Equillibration
equi_folder = os.path.join(DATA_OUT, 'sys_equi')
HA_time = 0.5
CA_time = 1.0
CA_LOW_time = 4.0

dt_HA = 0.001
dt = 0.002

HA_step = 1000 * HA_time / dt_HA
CA_step = 1000 * CA_time / dt
CA_LOW_step = 1000 * CA_LOW_time / dt

#Load from checkpoint
with open('checkpoint_em_20240915.pycpt', 'rb') as py_cpt:
    complex_sys = pickle.load(py_cpt)

#Parallelization
nthreads = int(os.environ.get('PBS_NCPUS', '12'))
#Set Parallelization
complex_sys.nt = nthreads
#complex_sys.ntmpi = 1
complex_sys.gpu_id = '0'

complex_sys.em_equi_three_step_iter_error(out_folder=equi_folder,
    no_constr_nsteps=em_step_number,
    constr_nsteps=em_step_number,
    nsteps_HA=HA_step,  
    nsteps_CA=CA_step,
    nsteps_CA_LOW=CA_LOW_step,
    dt=dt, dt_HA=dt_HA,
    tc_grps='Protein HZ1 Water_and_ions',
    tau_t= '0.1 0.1 0.1',
    ref_t= '310 310 310',
    vsite=vsite, maxwarn=10, iter_num=1)

#Checkpoint again
with open('checkpoint_equi_20240915.pycpt', 'wb') as py_cpt:
    pickle.dump(complex_sys, py_cpt)