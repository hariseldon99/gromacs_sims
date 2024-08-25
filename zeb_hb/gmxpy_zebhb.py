#!/usr/bin/env python
#In VSCODE, select python interpreter with ctrl+shift+p

from gromacs_py import gmx
import os

sys_name = "zeb_hb"
pdb = "zeb_hb/ZEB_HB_Refined.pdb"
DATA_OUT = "zeb_hb/tmp"

#Parallelization
nthreads = 4



#Creation of topology and posre files
md_sys = gmx.GmxSys(name=sys_name, coor_file=pdb)
#Set Parallelization and GPU
md_sys.nt = nthreads
md_sys.ntmpi = 1
#md_sys.ntomp = nthreads
md_sys.gpu_id = 0


#Hydrogen atoms need to be ignored, or else this won'AssertionErrort work with this particular pdb
md_sys.add_top(out_folder=DATA_OUT, name=sys_name, pdb2gmx_option_dict={'ignh': None})

md_sys.create_box(dist=1.0, box_type="dodecahedron", check_file_out=True)

#solvate and add ions
md_sys.solvate_add_ions(out_folder=DATA_OUT, name=sys_name, ion_C=0.15, box_dist=1.0)

md_sys.display()

#Energy Minimization
#Minimize a pdb structure in 2 steps, the first step without bonds constraints and the second step with
md_sys.em_2_steps(out_folder=DATA_OUT+"/em", name=sys_name, no_constr_nsteps=1000, constr_nsteps=1000, posres="", create_box_flag=False)
md_sys.convert_trj(traj=False)


md_sys.display()