#!/usr/bin/env python
import time
import os, sys
import glob
import subprocess
from biobb_gromacs.gromacs.mdrun import mdrun
#from biobb_gromacs.gromacs.grompp import grompp
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str

ext_time_ns = 40 # Extension time in nanoseconds
numprocs = 24 # Number of processors to use
mpithreads = 1 # Number of MPI threads to use
usegpu = True  # Whether to use GPU acceleration
gpuid = '0'  # GPU id to use, '0' for the first GPU, '1' for the second, etc.

simulation_dirnames = ["aminostep2012C",  
"aminotob2016C",  
"levparc2018C_chainB", 
"ndxgyra2018C",         
"ndxparc2018C_chainA", 
"pbppen1999C",
"aminostep2016C",
"levgyra2018C", 
"ndxgyra1996C", 
"ndxparc1999C_chainB", 
"pbpmez2018C",
"pbppen2018C"]

# ## Unset OMP_NUM_THREADS
    # 
    # The environment variable `OMP_num_threads_omp` is often set automatically by grid engines even if the user did not. Doing so creates conflicts with thread setting variables whenever `gromacs` runs. However, unsetting it before launching this notebook will allow `acpype` (the AMBER-GAFF2 topology generator) to use all available cores in the runtime machine, which might violate grid policy. So best to unset it after `acpype` runs.
    # Backup the value of OMP_NUM_THREADS to a string
omp_num_threads_backup = os.environ.get('OMP_NUM_THREADS', None)
if omp_num_threads_backup is not None:
    del os.environ['OMP_NUM_THREADS']

for simdir in simulation_dirnames:

    current_working_directory = os.getcwd()
    target_dir = os.path.join(current_working_directory, simdir)
    if not os.path.isdir(target_dir):
        raise FileNotFoundError(f"Directory '{simdir}' does not exist in {current_working_directory}")
        exit(1)
    os.chdir(target_dir)
    print(f"Entering directory {simdir}")

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    log_filename = f"extended_simulation_log_{timestamp}.log"
    time.sleep(5)
    sys.stdout = open(log_filename, "w")
    sys.stderr = sys.stdout

    print("=" * 50)
    print("Beginning of Extended Protein-Ligand Dynamics Simulation")
    print("=" * 50)

    #Assuming an nvt-npt equilibration, find the gro file from the last equilibration step
    input_npt_gro_files = glob.glob(f"*{simdir}*_npt.gro")
    # Find the full system topology file (after boxing, solvation, ion addition)
    genion_top_zip_files = glob.glob(f"*{simdir}*_genion_top.zip")
    #Now for the index file
    input_ndx_files = glob.glob(f"*{simdir}*_index.ndx")
    #Now, find the load the binary topology file
    input_tpr_files = glob.glob(f"*{simdir}*_gppmd.tpr")
    #Finally, get the checkpoint file
    input_cpt_files = glob.glob(f"*{simdir}*_md.cpt")

    #Now, for naming the extended sim outputs
    output_trr_files = glob.glob(f"*{simdir}*_md.trr")
    output_edr_files = glob.glob(f"*{simdir}*_md.edr")
    output_log_files = glob.glob(f"*{simdir}*_md.log")
    output_xtc_files = glob.glob(f"*{simdir}*_md.xtc")
    output_gro_files = glob.glob(f"*{simdir}*_md.gro")
    output_cpt_files = glob.glob(f"*{simdir}*_md.cpt")
    output_dry_gro_files = glob.glob(f"*{simdir}*_md_dry.gro")

    try:
        input_npt_gro = input_npt_gro_files[0]
        print(f"Found NPT equilibration gro file: {input_npt_gro}")
        output_gro = input_npt_gro.removesuffix('.gro')+'_ext.gro'
        print(f"Extended output gro file will be: {output_gro}")

        genion_top_zip = genion_top_zip_files[0]
        print(f"Found topology zip file: {genion_top_zip}")
        
        
        input_ndx = input_ndx_files[0]
        print(f"Found index file: {input_ndx}")

        input_tpr = input_tpr_files[0]
        print(f"Found input tpr file: {input_tpr}")

        input_cpt = input_cpt_files[0]
        print(f"Found checkpoint file: {input_cpt}")

        output_trr = output_trr_files[0]
        print(f"Found output trr file: {output_trr}")
        output_trr = output_trr.removesuffix('.trr')+'_ext.trr'
        print(f"Extended output trr file will be: {output_trr}")

        output_edr = output_edr_files[0]
        print(f"Found output edr file: {output_edr}")
        output_edr = output_edr.removesuffix('.edr')+'_ext.edr'
        print(f"Extended output edr file will be: {output_edr}")
        output_log = output_log_files[0]
        print(f"Found output log file: {output_log}")
        output_log = output_log.removesuffix('.log')+'_ext.log'
        print(f"Extended output log file will be: {output_log}")
        output_xtc = output_xtc_files[0]
        print(f"Found output xtc file: {output_xtc}")
        output_xtc = output_xtc.removesuffix('.xtc')+'_ext.xtc'
        print(f"Extended output xtc file will be: {output_xtc}")
        output_cpt = output_cpt_files[0]
        print(f"Found output cpt file: {output_cpt}")
        output_cpt = output_cpt.removesuffix('.cpt')+'_ext.cpt'
        print(f"Extended output cpt file will be: {output_cpt}")
        output_dry_gro = output_dry_gro_files[0]
        print(f"Found output dry gro file: {output_dry_gro}")
        output_dry_gro = output_dry_gro.removesuffix('_dry.gro')+'_ext_dry.gro'
        print(f"Extended output dry gro file will be: {output_dry_gro}")
    
    except FileNotFoundError as e:
        print(f"Error: {e}")
        exit(1)


    # -------------------------
    # Step 1: Extend simulation by 40 ns (40000 ps)
    # -------------------------
    # Use gmx convert-tpr (subprocess call), Biobb has no wrapper for extending
    extended_tpr = input_tpr.removesuffix('.tpr')+'_ext.tpr'
    subprocess.run([
        "gmx", "convert-tpr",
        "-s", input_tpr,
        "-extend", str(ext_time_ns * 1000),
        "-o", extended_tpr
    ], check=True)
    print(f"Extended TPR file created: {extended_tpr}")

    mdrun(input_tpr_path=extended_tpr,
        input_cpt_path=input_cpt,
        output_trr_path=output_trr,
        output_gro_path=output_gro,
        output_edr_path=output_edr,
        output_log_path=output_log,
        output_xtc_path=output_xtc,
        output_cpt_path=output_cpt,
        use_gpu=usegpu,
        num_threads_omp=numprocs,
        num_threads_mpi=mpithreads,
        gpu_id=gpuid)

    # Create prop dict and inputs/outputs (my way, which worked better)
    output_cluster_traj = output_xtc.removesuffix('.xtc')+'_cluster_traj.xtc'
    prop = {
        'cluster_selection':  'Protein_Other',
        'output_selection': 'Protein_Other',
        'center' : False,
        'center_selection': 'Protein_Other', #This is due to a bug in biobb_gromacs     
        'pbc' : 'cluster'
    }

    # Create and launch bb 
    gmx_image(input_traj_path=output_xtc,
            input_top_path=extended_tpr,
            input_index_path=input_ndx,
            output_traj_path=output_cluster_traj, 
            properties=prop)
    
    output_center_traj = output_xtc.removesuffix('.xtc')+'_cluster_center_traj.xtc'
    prop = {
        'center_selection':  'Protein_Other',
        'output_selection': 'Protein_Other',
        'pbc' : 'none',
        'center' : True
    }
    # Create and launch bb 
    gmx_image(input_traj_path=output_cluster_traj,
            input_top_path=extended_tpr,
            input_index_path=input_ndx,
            output_traj_path=output_center_traj, 
            properties=prop)


    # <a id="ppStep2"></a>
    # ### Step 2: Generating the output *dry* structure.
    # **Removing water molecules and ions** from the resulting structure

    # GMXTrjConvStr: Converting and/or manipulating a structure
    #                Removing water molecules and ions from the resulting structure
    #                The "dry" structure will be used as a topology to visualize 
    #                the "imaged dry" trajectory generated in the previous step.

    # Create prop dict and inputs/outputs
    prop = {
        'selection':  'Protein_Other'
    }

    # Create and launch bb
    gmx_trjconv_str(input_structure_path=output_gro,
                    input_top_path=extended_tpr,
                    input_index_path=input_ndx,
                    output_str_path=output_dry_gro, 
                    properties=prop)
    print(f"Extended MD simulation completed for {simdir}")
    os.chdir(current_working_directory)
    # Reset standard output and standard error back to default
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    print("=" * 50)
    print("End of Extended Protein-Ligand Dynamics Simulation")
    print("=" * 50)
    
#Before exiting, restore the environment variable OMP_NUM_THREADS
if omp_num_threads_backup is not None:
    os.environ['OMP_NUM_THREADS'] = omp_num_threads_backup