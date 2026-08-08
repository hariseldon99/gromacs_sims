#!/usr/bin/env python
"""
# Protein Ligand Complex MD Setup based on tutorial using BioExcel Building Blocks (biobb)
### Based on the official Gromacs tutorial: http://www.mdtutorials.com/gmx/complex/index.html
***
This tutorial aims to illustrate the process of **setting up a simulation system** containing a **protein in complex with a ligand**, step by step, using the **BioExcel Building Blocks library (biobb)**. 
***
**Biobb modules** used:

 - [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
 - [biobb_model](https://github.com/bioexcel/biobb_model): Tools to model macromolecular structures.
 - [biobb_chemistry](https://github.com/bioexcel/biobb_chemistry): Tools to manipulate chemical data.
 - [biobb_gromacs](https://github.com/bioexcel/biobb_gromacs): Tools to setup and run Molecular Dynamics simulations.
 - [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.
 - [biobb_structure_utils](https://github.com/bioexcel/biobb_structure_utils):  Tools to modify or extract information from a PDB structure file.
 
***
### Pipeline steps:
 16. [Post-processing Resulting 3D Trajectory](#post)
"""
import argparse
import os
import sys

from biobb_gromacs.gromacs.make_ndx import make_ndx

from biobb_analysis.gromacs.gmx_image import gmx_image

from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from tabulate import tabulate
import time


def molecular_dynamics_pp(complex):
    """    
    Main function to set up a molecular dynamics simulation PP for a protein-ligand complex dynamics.
    This function performs the following steps:
    1. Post-processes the resulting trajectory.
    """
    # Redirect standard output and error to a log file with a timestamped name
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    outdir = complex['outdir']
    os.makedirs(outdir, exist_ok=True)
    log_filename = os.path.join(complex['outdir'], f"simulation_PP_log_{timestamp}.log")
    
    # Tabulate the complex dictionary
    complex_table = [
        ["Output Directory", complex['outdir']],
        ["Log File", log_filename],
        ["gppmd_tpr", complex['gppmd_tpr']],
       # ["min_gro", complex['min_gro']],
        ["md_gro", complex['md_gro']],
        ["md_xtc", complex['md_xtc']],
        ["complex_ndx", complex['complex_ndx']]
    ]
    print("=" * 50)
    print("Start of Protein-Ligand Dynamics Post-Process")
    print("=" * 50)
    print(tabulate(complex_table, headers=["Parameter", "Value"], tablefmt="grid"))
    
    
    # Pause for 5 seconds
    time.sleep(5)
    sys.stdout = open(log_filename, "w")
    sys.stderr = sys.stdout
    
    # Store the current working directory
    current_working_directory = os.getcwd()
    os.chdir(outdir)
    
    # Create prop dict and inputs/outputs
    output_gppmd_tpr = complex['gppmd_tpr']
    #output_min_gro = complex['min_gro']
    output_md_gro = complex['md_gro']
    output_md_xtc = complex['md_xtc']
    output_complex_ndx = complex['complex_ndx']

    # Recreate selection indices
    # Create prop dict and inputs/outputs
    #prop = {
    #    'selection': "\"Protein\"|\"Other\"" 
    #}

    # Create and launch bb
    #make_ndx(input_structure_path=output_min_gro,
    #        output_ndx_path=output_complex_ndx,
    #        properties=prop)

    
    # Create prop dict and inputs/outputs (my way, which worked better)
    output_cluster_traj = complex['md_xtc'].removesuffix('.xtc')+'_cluster_traj.xtc'
    prop = {
        'cluster_selection':  'Protein_Other',
        'output_selection': 'Protein_Other',
        'center' : False,
        'center_selection': 'Protein_Other', #This is due to a bug in biobb_gromacs     
        'pbc' : 'cluster'
    }

    # Create and launch bb 
    gmx_image(input_traj_path=output_md_xtc,
            input_top_path=output_gppmd_tpr,
            input_index_path=output_complex_ndx,
            output_traj_path=output_cluster_traj, 
            properties=prop)

    output_center_traj = complex['md_xtc'].removesuffix('.xtc')+'_cluster_center_traj.xtc'
    prop = {
        'center_selection':  'Protein_Other',
        'output_selection': 'Protein_Other',
        'pbc' : 'none',
        'center' : True
    }
    # Create and launch bb 
    gmx_image(input_traj_path=output_cluster_traj,
            input_top_path=output_gppmd_tpr,
            input_index_path=output_complex_ndx,
            output_traj_path=output_center_traj, 
            properties=prop)


    # <a id="ppStep2"></a>
    # ### Step 2: Generating the output *dry* structure.
    # **Removing water molecules and ions** from the resulting structure

    output_dry_gro = complex['md_gro'].removesuffix('.gro')+'_md_dry.gro'
    prop = {
        'selection':  'Protein_Other'
    }

    # Create and launch bb
    gmx_trjconv_str(input_structure_path=output_md_gro,
                    input_top_path=output_gppmd_tpr,
                    input_index_path=output_complex_ndx,
                    output_str_path=output_dry_gro, 
                    properties=prop) 
    
    os.chdir(current_working_directory)
    # Reset standard output and standard error back to default
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    print("=" * 50)
    print("End of Protein-Ligand Dynamics PP")
    print("=" * 50)
    return True

if __name__ == '__main__':
    # Execute main function
    parser = argparse.ArgumentParser(description="Protein-Ligand Molecular Dynamics PP Setup")
    parser.add_argument("--outdir", type=str, default="./outputs", help="Output directory (default: ./outputs)")
    parser.add_argument("--gppmd_tpr", type=str, required=True, help="Path to the gppmd tpr file")
    #parser.add_argument("--min_gro", type=str, required=True, help="Path to the minimized gro file")
    parser.add_argument("--md_gro", type=str, required=True, help="Path to the MD gro file")
    parser.add_argument("--md_xtc", type=str, required=True, help="Path to the MD xtc file")
    parser.add_argument("--complex_ndx", type=str, required=True, help="Path to the selection index file")

    args = parser.parse_args()
    
    complex = {
        'outdir': args.outdir,
        'gppmd_tpr': args.gppmd_tpr,
        #"min_gro": args.min_gro,
        "md_gro": args.md_gro,
        "md_xtc": args.md_xtc,
        "complex_ndx": args.complex_ndx,
    }


    molecular_dynamics_pp(complex)